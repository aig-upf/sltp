#!/usr/bin/env python

#  Copyright (C) 2018-<date> Blai Bonet
#
#  Permission is hereby granted to distribute this software for
#  non-commercial research purposes, provided that this copyright
#  notice is included with any such distribution.
#
#  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
#  EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
#  SOFTWARE IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
#  ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#  Blai Bonet, bonet@ldc.usb.ve, bonetblai@gmail.com
#  Guillem Frances, guillem.frances@unibas.ch
import itertools
import logging
from enum import Enum

import numpy as np

from .util.tools import minimize_dnf
from .errors import CriticalPipelineError
from .util.console import header, lines
from .util.cnfwriter import CNFWriter
from .solvers import solve
from .util.performance import print_memory_usage
from .returncodes import ExitCode


class OptimizationPolicy(Enum):
    NUM_FEATURES = 1  # Minimize number of features
    TOTAL_FEATURE_COMPLEXITY = 2  # Minimize the sum of depths of selected features
    NUMERIC_FEATURES_FIRST = 3   # Minimize number of numeric features first, then overall features.


def compute_d2_index(s0, s1, t0, t1):
    assert s0 != t0
    if s0 < t0:
        return s0, s1, t0, t1
    else:
        return t0, t1, s0, s1


def generate_maxsat_problem(config, data, rng):
    translator, sample = create_maxsat_translator(config, data.sample)
    result = translator.run(config.cnf_filename, data.enforced_feature_idxs, data.in_goal_features)
    # translator.writer.print_variables(config.maxsat_variables_file)

    # TODO Serialize less stuff. Probably we don't need to serialize the full translator
    data = dict(cnf_translator=translator, post_cnf_sample=sample) if result == ExitCode.Success else dict()
    return result, data


class CompletenessInfo:
    def __init__(self, all_states, all_transitions, optimal_states, optimal_transitions):
        self.all_states = all_states
        self.all_transitions = all_transitions
        self.optimal_states = optimal_states
        self.optimal_transitions = optimal_transitions


def compute_completeness_info(sample, complete_only_wrt_optimal):
    """ Compute optimal states and transitions based on the experiment configuration """
    opt = complete_only_wrt_optimal
    all_states = set(sample.get_sorted_state_ids())
    all_transitions = [(s, s_prime) for s in sample.transitions for s_prime in sample.transitions[s]]

    if opt:
        optimal_states = sample.compute_optimal_states(include_goals=False)
        optimal_transitions = sample.optimal_transitions
    else:
        optimal_states = all_states.copy()
        optimal_transitions = set((x, y) for x, y in all_transitions)  # Consider all transitions as optimal

    cinfo = CompletenessInfo(all_states, all_transitions, optimal_states, optimal_transitions)
    return cinfo


def create_maxsat_translator(config, sample):
    optimization = config.optimization if hasattr(config, "optimization") else OptimizationPolicy.NUM_FEATURES

    feature_complexity = np.load(config.feature_complexity_filename + ".npy")
    feature_names = np.load(config.feature_names_filename + ".npy")
    feat_matrix = np.load(config.feature_matrix_filename + ".npy")
    bin_feat_matrix = np.load(config.bin_feature_matrix_filename + ".npy")
    feature_types = feat_matrix.max(0) <= 1  # i.e. feature_types[i] is True iff i is an (empirically) bool feature

    if not sample.goals:
        raise CriticalPipelineError("No goal state identified in the sample, SAT theory will be trivially satisfiable")

    cinfo = compute_completeness_info(sample, config.complete_only_wrt_optimal)

    # Remove states that are redundant
    if config.prune_redundant_states:
        sample = preprocess_sample(sample, feat_matrix, bin_feat_matrix, cinfo)

        # Reproject the denotation matrices to the new state indices
        assert sample.remapping
        projection = list(sorted(sample.remapping.keys()))
        feat_matrix = feat_matrix[projection]
        bin_feat_matrix = bin_feat_matrix[projection]
        cinfo = compute_completeness_info(sample, config.complete_only_wrt_optimal)
    else:
        logging.warning("Not pruning redundant states")

    translator = ModelTranslator(feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types,
                                 sample, optimization, cinfo)
    return translator, sample


def run_solver(config, data, rng):
    solution = solve(config.experiment_dir, config.cnf_filename, config.maxsat_solver, config.maxsat_timeout)
    if not solution.solved and solution.result == "UNSATISFIABLE":
        return ExitCode.MaxsatModelUnsat, dict()
    else:
        logging.info("MAXSAT solution with cost {} found".format(solution.cost))

    return ExitCode.Success, dict(cnf_solution=solution)


class ModelTranslator:
    def __init__(self, feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types,
                 sample, optimization, cinfo):

        self.writer = None
        self.feat_matrix = feat_matrix
        self.bin_feat_matrix = bin_feat_matrix
        self.feature_complexity = feature_complexity
        self.feature_names = feature_names
        self.feature_types = feature_types
        self.sample = sample
        self.cinfo = cinfo

        self.non_goal_states = cinfo.all_states.difference(self.sample.goals)
        self.d1_pairs = self.compute_d1_pairs()

        # Compute for each pair s and t of states which features distinguish s and t
        self.d1_distinguishing_features = self.compute_d1_distinguishing_features(bin_feat_matrix)

        self.var_selected = None
        self.var_d2 = None

        self.n_selected_clauses = 0
        self.n_d1_clauses = 0
        self.n_d2_clauses = 0
        self.n_bridge_clauses = 0
        self.n_goal_clauses = 0
        self.qchanges = dict()
        self.compute_feature_weight = self.setup_optimization_policy(optimization)

    def compute_d1_pairs(self):
        pairs = set(self.iterate_over_bridge_pairs())  # initialize with all bridge pairs
        for (x, y) in itertools.product(self.sample.goals, self.cinfo.all_states):
            if y in self.sample.goals and y <= x:
                # Break symmetries: if both x and y are goals, only need to take into account 1 of the 2 permutations
                continue
            pairs.add((x, y))
        return pairs
        # return list(x for x in pairs)

    def d1idx(self, s, t):
        assert (s, t) in self.d1_pairs or (t, s) in self.d1_pairs
        return (s, t) if (s, t) in self.d1_pairs else (t, s)

    def iterate_over_bridge_pairs(self):
        for (x, y) in itertools.product(self.cinfo.optimal_states, self.cinfo.all_states):
            if y in self.cinfo.optimal_states and y <= x:
                # Break symmetries: if both x and y optimal, only need to take into account 1 of the 2 permutations
                continue
            yield (x, y)

    def compute_d1_distinguishing_features(self, bin_feat_matrix):
        """ """
        nf = bin_feat_matrix.shape[1]
        npairs = len(self.d1_pairs)
        nentries = nf * npairs  # How many feature denotations we'll have to compute

        logging.info("Computing sets of D1-distinguishing features for {} state pairs "
                     "and {} features ({:0.1f}M matrix entries)".format(npairs, nf, nentries / 1000000))
        # We extract the non-zero indexes of the XOR of both binary matrices and store them in a set:
        return {(s1, s2): set(np.nonzero(np.logical_xor(bin_feat_matrix[s1], bin_feat_matrix[s2]))[0].flat)
                for s1, s2 in self.d1_pairs}

    def create_bridge_clauses(self, s_t_distinguishing, s, t):
        # If there are no transitions (t, t') in the sample set, then we do not post the bridge constraint.
        # if t not in self.transitions or not self.transitions[t]:
        # TODO WE MIGHT NEED A BETTER WAY OF IDENTIFYING NON-EXPANDED STATES AND DISTINGUISHING THEM FROM
        # TODO STATES WHICH HAVE NO SUCCESSOR
        if not self.sample.transitions[s] or not self.sample.transitions[t]:
            return

        assert s != t
        assert s in self.cinfo.optimal_states

        for s_prime in self.sample.transitions[s]:
            if not self.is_optimal_transition(s, s_prime):
                continue

            # print("Bridge clause btw transitions {} and {}".format((s, s_prime), (t, "?")))
            # Start with OR_i Selected(f_i)
            lits = [self.writer.literal(self.var_selected[f], True) for f in s_t_distinguishing]
            for t_prime in self.sample.transitions[t]:
                idx = compute_d2_index(s, s_prime, t, t_prime)
                # And add a literal -D2(s,s',t,t') for each child t' of t
                lits.append(self.writer.literal(self.var_d2[idx], False))
            self.writer.clause(lits)
            self.n_bridge_clauses += 1

    def retrieve_possibly_cached_qchanges(self, s0, s1):
        try:
            return self.qchanges[(s0, s1)]
        except KeyError:
            self.qchanges[(s0, s1)] = val = compute_qualitative_changes(self.feat_matrix, s0, s1)
            return val

    def is_optimal_transition(self, s, sprime):
        return (s, sprime) in self.cinfo.optimal_transitions

    def iterate_over_optimal_transitions(self):
        """ """
        return self.cinfo.optimal_transitions

    def iterate_over_all_transitions(self):
        for s in self.sample.transitions:
            for s_prime in self.sample.transitions[s]:
                yield (s, s_prime)

    def iterate_over_d2_transitions(self):
        for s, s_prime in self.iterate_over_optimal_transitions():
            for t, t_prime in self.iterate_over_all_transitions():
                if s == t or (self.is_optimal_transition(t, t_prime) and t < s):  # Symmetry-breaking
                    continue
                yield (s, s_prime, t, t_prime)

    def create_variables(self):
        self.var_selected = \
            [self.writer.variable("F{}".format(feat)) for feat in range(0, self.feat_matrix.shape[1])]

        self.var_d2 = dict()

        for s, s_prime, t, t_prime in self.iterate_over_d2_transitions():
            idx = compute_d2_index(s, s_prime, t, t_prime)
            varname = "D2[({},{}), ({},{})]".format(*idx)
            self.var_d2[idx] = self.writer.variable(varname)

        self.writer.close()
        logging.info("A total of {} MAXSAT variables were created".format(len(self.writer.variables)))

    def run(self, cnf_filename, enforced_feature_idxs, in_goal_features):
        logging.info("Generating MAXSAT problem from {} states, {} transitions and {} features"
                     .format(self.feat_matrix.shape[0], self.sample.num_transitions(), self.feat_matrix.shape[1]))

        # allf = list(range(0, len(self.feature_names)))  # i.e. select all features
        # in_goal_features = set()
        self.writer = CNFWriter(cnf_filename)
        self.create_variables()

        logging.info("Generating bridge constraints")
        for s1, s2 in self.iterate_over_bridge_pairs():
            # Only if D1(s, t) _might_ be false, we create the bridge clauses between values of D1 and D2:
            # For each s' in succ(s), we'll post:
            #    -D1(s,t) --> OR_t'  D2(s,s',t,t')
            # where t' ranges over all successors of state t
            if sum(1 for x in (s1, s2) if x in self.sample.goals) != 1:
                assert s1 in self.cinfo.optimal_states
                s_t_distinguishing = self.d1_distinguishing_features[self.d1idx(s1, s2)]
                self.create_bridge_clauses(s_t_distinguishing, s1, s2)
                if s2 in self.cinfo.optimal_states:
                    self.create_bridge_clauses(s_t_distinguishing, s2, s1)

        # Force D1(s1, s2) to be true if exactly one of the two states is a goal state
        logging.info("Generating goal constraints for {} state pairs".format(
            len(self.sample.goals)*len(self.non_goal_states)))
        for s1, s2 in itertools.product(self.sample.goals, self.non_goal_states):

            if len(in_goal_features):
                # we are enforcing goal-distinguishing features elsewhere, no need to do anything here
                if not in_goal_features.intersection(self.d1_distinguishing_features[self.d1idx(s1, s2)]):
                    raise RuntimeError(undist_goal_warning(s1, s2))

            else:
                if len(self.d1_distinguishing_features[self.d1idx(s1, s2)]) == 0:
                    logging.warning(undist_goal_warning(s1, s2))
                    return ExitCode.MaxsatModelUnsat

                # print("Force goal distinguishability: {}".format((s1, s2)))
                self.writer.clause([self.writer.literal(self.var_selected[f], True)
                                    for f in self.d1_distinguishing_features[self.d1idx(s1, s2)]])
                self.n_goal_clauses += 1

        logging.info("Generating D2 constraints from {} D2 variables".format(len(self.var_d2)))
        for s0, s1, t0, t1 in self.iterate_over_d2_transitions():
            d2_var = self.var_d2[compute_d2_index(s0, s1, t0, t1)]

            qchanges_s0s1 = self.retrieve_possibly_cached_qchanges(s0, s1)
            qchanges_t0t1 = self.retrieve_possibly_cached_qchanges(t0, t1)
            equal_idxs = np.equal(qchanges_s0s1, qchanges_t0t1)
            np_d2_distinguishing_features = np.where(equal_idxs == False)[0]

            np_d1_dist = self.d1_distinguishing_features[self.d1idx(s0, t0)]
            # D2(s0,s1,t0,t2) <-- OR_f selected(f), where f ranges over features that d2-distinguish the transition
            # but do _not_ d1-distinguish the two states at the origin of each transition.
            d2_lit = self.writer.literal(d2_var, True)
            # forward_clause_literals = [self.writer.literal(d2_var, False)]
            for f in np_d2_distinguishing_features.flat:
                if f not in np_d1_dist:
                    # forward_clause_literals.append(self.writer.literal(self.var_selected[f], True))
                    self.writer.clause([d2_lit, self.writer.literal(self.var_selected[f], False)])
                    self.n_d2_clauses += 1

            # self.writer.clause(forward_clause_literals)
            # self.n_d2_clauses += 1

        # Add the weighted clauses to minimize the number of selected features
        for feature, var in enumerate(self.var_selected):
            d = self.compute_feature_weight(feature)
            self.writer.clause([self.writer.literal(var, False)], weight=d)
            self.n_selected_clauses += 1

        # Make sure those features we want to enforce in the solution, if any, are posted as necessarily selected
        # for enforced in enforced_feature_idxs:
        #     self.writer.clause([self.writer.literal(self.var_selected[enforced], True)])
        assert not enforced_feature_idxs

        logging.info("Enforcing features: {}".format(", ".join(str(self.var_selected[x]) for x in in_goal_features)))
        for enforced in in_goal_features:
            self.writer.clause([self.writer.literal(self.var_selected[enforced], True)])

        # from pympler.asizeof import asizeof
        # print(f"asizeof(self.writer): {asizeof(self.writer)/(1024*1024)}MB")
        # print(f"asizeof(self): {asizeof(self)/(1024*1024)}MB")
        # print(f"asizeof(self.d1_distinguishing_features): {asizeof(self.d1_distinguishing_features)/(1024*1024)}MB")
        # print(f"np_d2_distinguishing_features.nbytes: {np_d2_distinguishing_features.nbytes/(1024*1024)}MB")
        print_memory_usage()

        self.report_stats()

        # self.debug_tests()
        self.d1_distinguishing_features = None  # We won't need this anymore
        self.writer.save()

        return ExitCode.Success

    def ftype(self, f):
        return bool if self.feature_types[f] else int

    def decode_solution(self, assignment):
        varmapping = self.writer.mapping
        true_variables = set(varmapping[idx] for idx, value in assignment.items() if value is True)
        feature_mapping = {variable: feature for feature, variable in enumerate(self.var_selected)}
        assert len(feature_mapping) == len(self.var_selected)
        selected_features = sorted(feature_mapping[v] for v in true_variables if v in feature_mapping)
        return selected_features

    def report_stats(self):

        d1_distinguishing = self.d1_distinguishing_features.values()
        avg_num_d1_dist_features = sum(len(feats) for feats in d1_distinguishing) / len(d1_distinguishing)
        n_undistinguishable_state_pairs = sum(1 for x in d1_distinguishing if len(x) == 0)
        s = header("Max-sat encoding stats", 1)
        s += "\n" + lines("Number of D1-undistinguishable state pairs: {}".format(n_undistinguishable_state_pairs), 1)
        s += "\n" + lines("Avg. # of D1-distinguishing features: {:.2f}".format(avg_num_d1_dist_features), 1)
        s += "\n" + lines("Clauses (possibly with repetitions):".format(), 1)
        s += "\n" + lines("Selected: {}".format(self.n_selected_clauses), 2)
        s += "\n" + lines("D1: {}".format(self.n_d1_clauses), 2)
        s += "\n" + lines("D2: {}".format(self.n_d2_clauses), 2)
        s += "\n" + lines("Bridge: {}".format(self.n_bridge_clauses), 2)
        s += "\n" + lines("Goal: {}".format(self.n_goal_clauses), 2)
        s += "\n" + lines("TOTAL: {}".format(self.n_selected_clauses + self.n_d1_clauses + self.n_d2_clauses +
                                             self.n_bridge_clauses + self.n_goal_clauses), 2)
        logging.info("\n{}".format(s))

    def debug_tests(self):
        # undist = self.undistinguishable_state_pairs
        # for s, t in undist:
        #     for s_prime in self.transitions[s]:
        #         some_d2_indistinguishable = False
        #         for t_prime in self.transitions[t]:
        #             idx = compute_d2_index(s, s_prime, t, t_prime)
        #             if idx not in self.d2_distinguished:
        #                 some_d2_indistinguishable = True
        #
        #         if not some_d2_indistinguishable:
        #             print("State pair {} violates bridge clauses".format((s, t)))
        #             for t_prime in self.transitions[t]:
        #                 for f in self.features:
        #                     print("Delta[{}] for feature {}: {}".format((s, s_prime),
        #                                                                 f, self.qchanges[(s, s_prime, f)]))
        #                     print("Delta[{}] for feature {}: {}".format((t, t_prime), f,
        #                                                                 self.qchanges[(t, t_prime, f)]))
        #             assert False
        pass

    def opt_policy_num_features(self, feature):
        return 1

    def opt_policy_total_feature_depth(self, feature):
        # return feature.complexity()
        return self.feature_complexity[feature]

    def setup_optimization_policy(self, optimization):
        if optimization == OptimizationPolicy.NUM_FEATURES:
            return self.opt_policy_num_features

        if optimization == OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY:
            return self.opt_policy_total_feature_depth

        raise RuntimeError("Unknown optimization policy")

    def print_qnp_precondition_atom(self, feature, value, namer):
        assert value in (True, False)
        return "{} {}".format(namer(self.feature_names[feature]), int(value))

    def print_qnp_conjunction(self, conjunction, namer):
        sorted_ = sorted(self.print_qnp_precondition_atom(f, v, namer) for f, v in conjunction)
        return "{} {}".format(len(sorted_), " ".join(sorted_))

    def compute_qnp(self, actions, features, config, sample):
        logging.info("Writing QNP abstraction to {}".format(config.qnp_abstraction_filename))
        namer = lambda x: config.feature_namer(x).replace(" ", "")  # let's make sure no spaces in feature names

        init_dnf, goal_dnf = self.compute_init_goals(sample, features)

        # TODO Just pick an arbitrary clause from the DNF. Should do better!
        init = next(iter(init_dnf))
        goal = next(iter(goal_dnf))

        if len(goal_dnf) > 1:
            logging.warning("Cannot minimize abstract goal DNF to single conjunction:")
            # _ = [print("\t" + self.print_qnp_conjunction(x, namer)) for x in goal_dnf]
        else:
            logging.info("Abstract goal is: {}".format(self.print_qnp_conjunction(goal, namer)))

        if len(init_dnf) > 1:
            logging.warning("Cannot minimize abstract initial state DNF to single conjunction:")
            # _ = [print("\t" + self.print_qnp_conjunction(x, namer)) for x in init_dnf]

        with open(config.qnp_abstraction_filename, "w") as f:
            print(config.domain_dir, file=f)

            name_types = " ".join("{} {}".format(namer(f["name"]), f["type"]) for f in features)
            feat_line = "{} {}".format(len(features), name_types)
            print(feat_line, file=f)

            print(self.print_qnp_conjunction(init, namer), file=f)
            print(self.print_qnp_conjunction(goal, namer), file=f)

            # Actions
            print("{}".format(len(actions)), file=f)
            for i, action in enumerate(actions, 1):
                print("action_{}".format(i), file=f)
                print(self.print_qnp_conjunction(action.preconditions, namer), file=f)

                effs = sorted(eff.print_qnp_named(namer) for eff in action.effects)
                print("{} {}".format(len(effs), " ".join(effs)), file=f)

    def extract_dnf(self, feature_indexes, states):
        """ """
        dnf = self.bin_feat_matrix[np.ix_(states, feature_indexes)]  # see https://stackoverflow.com/q/22927181
        dnf = set(frozenset((i, val) for i, val in zip(feature_indexes, f)) for f in dnf)
        return minimize_dnf(dnf)

    def compute_init_goals(self, sample, features):
        fidxs = [f["idx"] for f in features]

        init_dnf = self.extract_dnf(fidxs, sorted(sample.roots))
        goal_dnf = self.extract_dnf(fidxs, sorted(sample.goals))

        return init_dnf, goal_dnf


def compute_qualitative_changes(feature_matrix, s0, s1):
    return np.sign(feature_matrix[s1] - feature_matrix[s0])


def preprocess_sample(sample, feat_matrix, bin_feat_matrix, cinfo):
    logging.info("Preprocessing sample {} to prune redundant states".format(sample))

    nonoptimal_states = cinfo.all_states - cinfo.optimal_states
    isomorphisms = dict()

    for s, t in itertools.product(cinfo.optimal_states, cinfo.all_states):
        if s != t:
            check_isomorphic(sample, bin_feat_matrix, feat_matrix, s, t, isomorphisms)

    for s, t in itertools.combinations(nonoptimal_states, 2):
        check_isomorphic(sample, bin_feat_matrix, feat_matrix, s, t, isomorphisms)

    # logging.info("Set of prunable states ({}): {}".format(len(prunable_states), prunable_states))
    logging.info("Number of prunable states: {}".format(len(isomorphisms)))
    selected = cinfo.all_states - set(isomorphisms.keys())
    resampled = sample.resample(selected)
    logging.info("Processed sample: {}".format(resampled))
    return resampled


def has_analog_transition(sample, feat_matrix, s, t):
    """ Check whether all transitions starting in s have some transition starting in t with same qualitative nature
        on the set of all features in the given feature matrix
    """
    for sp in sample.transitions[s]:
        for tp in sample.transitions[t]:
            qchanges_s0s1 = compute_qualitative_changes(feat_matrix, s, sp)
            qchanges_t0t1 = compute_qualitative_changes(feat_matrix, t, tp)
            equal_idxs = np.equal(qchanges_s0s1, qchanges_t0t1)
            d2_distinguishing_features = np.where(equal_idxs == False)[0]
            if d2_distinguishing_features.size == 0:
                return True
    return False


def check_isomorphic(sample, bin_feat_matrix, feat_matrix, s, t, isomorphisms):
    if s in isomorphisms or t in isomorphisms:
        return
    distinguishing = np.nonzero(np.logical_xor(bin_feat_matrix[s], bin_feat_matrix[t]))[0]
    if distinguishing.size == 0:  # s and t are not distinguishable
        if sum(1 for x in (s, t) if x in sample.goals) == 1:  # Only one of the two is a goal!
            raise RuntimeError(undist_goal_warning(s, t))

        if has_analog_transition(sample, feat_matrix, s, t) and has_analog_transition(sample, feat_matrix, t, s):
            isomorphisms[t] = s
    return


def undist_goal_warning(s1, s2):
    return "No feature in pool can distinguish states {}, but one of them is a goal and the other is not." \
           " MAXSAT encoding will be UNSAT".format((s1, s2))
