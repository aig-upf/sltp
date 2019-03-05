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
import os
from collections import defaultdict
from enum import Enum

import numpy as np

from .sampling import TransitionSample
from .errors import CriticalPipelineError
from tarski.dl import FeatureValueChange, EmpiricalBinaryConcept
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
    translator = create_maxsat_translator(config, data)
    result = translator.run(config.cnf_filename, data.enforced_feature_idxs, data.in_goal_features)
    # translator.writer.print_variables(config.maxsat_variables_file)

    # TODO Serialize less stuff. Probably we don't need to serialize the full translator
    data = dict(cnf_translator=translator) if result == ExitCode.Success else dict()
    return result, data


def create_maxsat_translator(config, data):
    optimization = config.optimization if hasattr(config, "optimization") else OptimizationPolicy.NUM_FEATURES
    sample = data.sample

    feature_complexity = np.load(config.feature_complexity_filename + ".npy")
    feature_names = np.load(config.feature_names_filename + ".npy")
    feat_matrix = np.load(config.feature_matrix_filename + ".npy")
    bin_feat_matrix = np.load(config.bin_feature_matrix_filename + ".npy")
    feature_types = feat_matrix.max(0) <= 1  # i.e. feature_types[i] is True iff i is an (empirically) bool feature

    logging.info("Generating MAXSAT problem from {} states, {} transitions and {} features"
                 .format(feat_matrix.shape[0], sample.num_transitions(), feat_matrix.shape[1]))

    if not sample.goals:
        raise CriticalPipelineError("No goal state identified in the sample, SAT theory will be trivially satisfiable")

    translator = ModelTranslator(feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types,
                                 sample, optimization, config.complete_only_wrt_optimal)
    return translator


def run_solver(config, data, rng):
    solution = solve(config.experiment_dir, config.cnf_filename, config.maxsat_solver, config.maxsat_timeout)
    if not solution.solved and solution.result == "UNSATISFIABLE":
        return ExitCode.MaxsatModelUnsat, dict()
    else:
        logging.info("MAXSAT solution with cost {} found".format(solution.cost))

    return ExitCode.Success, dict(cnf_solution=solution)


class AbstractAction:
    def __init__(self, preconditions, effects, name=None):
        self.name = name
        self.preconditions = preconditions
        self.effects = frozenset(effects)
        self.hash = hash((self.__class__, self.preconditions, self.effects))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.preconditions == other.preconditions and self.effects == other.effects)

    def __str__(self):
        precs = " and ".join(map(str, self.preconditions))
        effs = ", ".join(map(str, self.effects))
        return "AbstractAction<{}; {}>".format(precs, effs)

    __repr__ = __str__


class ActionEffect:
    def __init__(self, feature, feature_name, change):
        self.feature = feature
        self.feature_name = feature_name
        self.change = change
        self.hash = hash((self.__class__, self.feature, self.change))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.feature == other.feature and self.change == other.change)

    def __str__(self):
        return self.print_named()

    __repr__ = __str__

    def print_named(self, namer=lambda s: s):
        name = namer(self.feature_name)
        if self.change == FeatureValueChange.ADD:
            return name
        if self.change == FeatureValueChange.DEL:
            return "NOT {}".format(name)
        if self.change == FeatureValueChange.ADD_OR_NIL:
            return "ADD* {}".format(name)
        if self.change == FeatureValueChange.INC:
            return "INC {}".format(name)
        if self.change == FeatureValueChange.DEC:
            return "DEC {}".format(name)
        if self.change == FeatureValueChange.INC_OR_NIL:
            return "INC* {}".format(name)
        raise RuntimeError("Unexpected effect type")

    def print_qnp_named(self, namer=lambda s: s):
        name = namer(self.feature_name)
        if self.change in (FeatureValueChange.ADD, FeatureValueChange.INC):
            return "{} 1".format(name)

        if self.change in (FeatureValueChange.DEL, FeatureValueChange.DEC):
            return "{} 0".format(name)

        if self.change in (FeatureValueChange.INC_OR_NIL, FeatureValueChange.ADD_OR_NIL):
            assert False, "Relaxed INC semantics not supported for QNP"
        raise RuntimeError("Unexpected effect type")


def prettyprint_atom(feature, polarity, namer):
    named = namer(feature)
    if isinstance(feature, EmpiricalBinaryConcept):
        return named if polarity else "not {}".format(named)
    return "{} > 0".format(named) if polarity else "{} = 0".format(named)


def prettyprint_abstract_action(action, all_features, namer):
    precs = " and ".join(map(lambda p: prettyprint_atom(all_features[p[0]], p[1], namer), action.preconditions))
    effs = ", ".join(map(lambda eff: eff.print_named(namer), action.effects))
    return "AbstractAction<{}; {}>".format(precs, effs)


class ModelTranslator:
    def __init__(self, feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types,
                 sample: TransitionSample,
                 optimization, complete_only_wrt_optimal):

        self.writer = None
        self.feat_matrix = feat_matrix
        self.bin_feat_matrix = bin_feat_matrix
        self.feature_complexity = feature_complexity
        self.feature_names = feature_names
        self.feature_types = feature_types
        self.complete_only_wrt_optimal = complete_only_wrt_optimal
        self.sample = sample
        self.state_ids = sample.get_sorted_state_ids()

        # i.e. all states:
        self.optimal_and_nonoptimal = set(sample.get_sorted_state_ids())

        self.optimal_transitions = sample.optimal_transitions

        if complete_only_wrt_optimal:
            self.optimal_states = sample.compute_optimal_states(include_goals=False)

        else:
            # Consider all transitions as optimal
            self.optimal_states = self.optimal_and_nonoptimal.copy()
            self.optimal_transitions = set((x, y) for x, y in self.iterate_over_all_transitions())

        self.non_goal_states = self.optimal_and_nonoptimal.difference(self.sample.goals)
        self.d1_pairs = self.compute_d1_pairs()

        # Compute for each pair s and t of states which features distinguish s and t
        self.d1_distinguishing_features = self.compute_d1_distinguishing_features(bin_feat_matrix)

        self.var_selected = None
        self.var_d1 = None
        self.var_d2 = None

        self.n_selected_clauses = 0
        self.n_d1_clauses = 0
        self.n_d2_clauses = 0
        self.n_bridge_clauses = 0
        self.n_goal_clauses = 0
        self.qchanges = dict()
        self.compute_feature_weight = self.setup_optimization_policy(optimization)

    def compute_d1_pairs(self):
        pairs = set()
        for (x, y) in itertools.product(self.optimal_states, self.optimal_and_nonoptimal):
            if y in self.optimal_states and y <= x:
                # Break symmetries: if both x and y optimal, only need to take into account 1 of the 2 permutations
                continue
            pairs.add((x, y))
        for (x, y) in itertools.product(self.sample.goals, self.non_goal_states):
            pairs.add((x, y))
        return list(x for x in pairs)

    def iterate_over_bridge_pairs(self):
        for (x, y) in itertools.product(self.optimal_states, self.optimal_and_nonoptimal):
            if y in self.optimal_states and y <= x:
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
        assert s in self.optimal_states

        for s_prime in self.sample.transitions[s]:
            if not self.is_optimal_transition(s, s_prime):
                continue

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
            self.qchanges[(s0, s1)] = val = self.compute_qchanges(s0, s1)
            return val

    def is_optimal_transition(self, s, sprime):
        return (s, sprime) in self.optimal_transitions

    def iterate_over_optimal_transitions(self):
        """ """
        return self.optimal_transitions

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
        # self.var_selected = {feat: self.writer.variable("selected({})".format(feat)) for feat in self.features}
        self.var_selected = \
            [self.writer.variable("selected(F{})".format(feat)) for feat in range(0, self.feat_matrix.shape[1])]

        self.var_d2 = dict()

        for s, s_prime, t, t_prime in self.iterate_over_d2_transitions():
            idx = compute_d2_index(s, s_prime, t, t_prime)
            varname = "D2[({},{}), ({},{})]".format(*idx)
            self.var_d2[idx] = self.writer.variable(varname)

        self.writer.close()
        logging.info("A total of {} MAXSAT variables were created".format(len(self.writer.variables)))

    def compute_d1_idx(self, s, t):
        assert s in self.optimal_states

        if t in self.optimal_states:  # Break symmetries when both states are optimal
            return min(s, t), max(s, t)
        return s, t

    def run(self, cnf_filename, enforced_feature_idxs, in_goal_features):
        # allf = list(range(0, len(self.feature_names)))  # i.e. select all features
        # self.find_inconsistency(allf)
        # in_goal_features = set()
        self.writer = CNFWriter(cnf_filename)
        self.create_variables()

        logging.info("Generating bridge constraints from {} D1 pairs".format(len(self.d1_pairs)))
        for s1, s2 in self.iterate_over_bridge_pairs():
            # Only if D1(s, t) _might_ be false, we create the bridge clauses between values of D1 and D2:
            # For each s' in succ(s), we'll post:
            #    -D1(s,t) --> OR_t'  D2(s,s',t,t')
            # where t' ranges over all successors of state t
            if sum(1 for x in (s1, s2) if x in self.sample.goals) != 1:
                assert s1 in self.optimal_states
                s_t_distinguishing = self.d1_distinguishing_features[(s1, s2)]
                self.create_bridge_clauses(s_t_distinguishing, s1, s2)
                if s2 in self.optimal_states:
                    self.create_bridge_clauses(s_t_distinguishing, s2, s1)

        # Force D1(s1, s2) to be true if exactly one of the two states is a goal state
        logging.info("Generating goal constraints for {} state pairs".format(
            len(self.sample.goals)*len(self.non_goal_states)))
        for s1, s2 in itertools.product(self.sample.goals, self.non_goal_states):

            if len(in_goal_features):
                # we are enforcing goal-distinguishing features elsewhere, no need to do anything here
                if not in_goal_features.intersection(self.d1_distinguishing_features[(s1, s2)]):
                    raise RuntimeError("No feature in pool can distinguish states {}, "
                                       "but forced features should be enough to distinguish them".format((s1, s2)))

            else:
                if len(self.d1_distinguishing_features[(s1, s2)]) == 0:
                    logging.warning("No feature in pool can distinguish states {}, but one of them is a goal and "
                                    "the other is not. MAXSAT encoding will be UNSAT".format((s1, s2)))
                    return ExitCode.MaxsatModelUnsat

                self.writer.clause([self.writer.literal(self.var_selected[f], True)
                                    for f in self.d1_distinguishing_features[(s1, s2)]])
                self.n_goal_clauses += 1

        logging.info("Generating D2 constraints from {} D2 variables".format(len(self.var_d2)))
        for s0, s1, t0, t1 in self.iterate_over_d2_transitions():
            d2_var = self.var_d2[compute_d2_index(s0, s1, t0, t1)]

            qchanges_s0s1 = self.retrieve_possibly_cached_qchanges(s0, s1)
            qchanges_t0t1 = self.retrieve_possibly_cached_qchanges(t0, t1)
            equal_idxs = np.equal(qchanges_s0s1, qchanges_t0t1)
            np_d2_distinguishing_features = np.where(equal_idxs == False)[0]

            np_d1_dist = self.d1_distinguishing_features[self.compute_d1_idx(s0, t0)]
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

    def compute_qchanges(self, s0, s1):
        return np.sign(self.feat_matrix[s1] - self.feat_matrix[s0])

    def ftype(self, f):
        return bool if self.feature_types[f] else int

    def generate_eff_change(self, feat, qchange):

        if self.ftype(feat) == bool:
            assert qchange in (-1, 1)
            return {-1: FeatureValueChange.DEL, 1: FeatureValueChange.ADD, 2: FeatureValueChange.ADD_OR_NIL}[qchange]

        # else type must be int
        assert qchange in (-1, 1, 2)
        return {-1: FeatureValueChange.DEC, 1: FeatureValueChange.INC, 2: FeatureValueChange.INC_OR_NIL}[qchange]

    def compute_state_abstractions(self, features):
        abstract_states = set()
        state_abstraction = dict()
        for state in self.state_ids:
            abstract = self.compute_abstract_state(features, state)
            abstract_states.add(abstract)
            state_abstraction[state] = abstract
        return abstract_states, state_abstraction

    def compute_abstraction(self, selected_features, namer):
        if not selected_features:
            raise CriticalPipelineError("Zero-cost maxsat solution - "
                                        "no action model possible, the encoding has likely some error")
        print("Features (total complexity: {}): ".format(
            sum(self.feature_complexity[f] for f in selected_features)))
        print('\t' + '\n\t'.join("{}. {} [k={}, id={}]".format(
            i, namer(self.feature_names[f]), self.feature_complexity[f], f) for i, f in enumerate(selected_features, 0)))

        # logging.debug("Features:\n{}".format(serialize_to_string([features[i] for i in selected_features])))
        selected_data = [dict(idx=f,
                              name=namer(self.feature_names[f]),
                              type=int(not self.feature_types[f]))
                         for f in selected_features]

        abstract_states, state_abstraction = self.compute_state_abstractions(selected_features)

        abstract_actions = set()
        for s, sprime in self.iterate_over_optimal_transitions():
            abstract_s = state_abstraction[s]
            abstract_sprime = state_abstraction[sprime]

            qchanges = self.retrieve_possibly_cached_qchanges(s, sprime)
            selected_qchanges = qchanges[selected_features]
            abstract_effects = [ActionEffect(i, self.feature_names[f], self.generate_eff_change(f, c))
                                # No need to record "NIL" changes:
                                for i, (f, c) in enumerate(zip(selected_features, selected_qchanges)) if c != 0]

            precondition_bitmap = frozenset(enumerate(abstract_s))
            abstract_actions.add(AbstractAction(precondition_bitmap, abstract_effects))
            if len(abstract_effects) == 0:
                msg = "Abstract no-op necessary [concrete: ({}, {}), abstract: ({}, {})]".format(
                    s, sprime, abstract_s, abstract_sprime)
                logging.warning(msg)

        logging.info("Abstract state space: {} states and {} actions".
                     format(len(abstract_states), len(abstract_actions)))
        return abstract_states, abstract_actions, selected_data

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

    def compute_abstract_state(self, features, state_id):
        bin_values = self.bin_feat_matrix[state_id, features]
        return tuple(map(bool, bin_values))

    def compute_action_model(self, selected_features, config):
        states, actions, selected_data = self.compute_abstraction(selected_features, config.feature_namer)
        self.print_actions(actions, os.path.join(config.experiment_dir, 'actions.txt'), config.feature_namer)
        states, actions = optimize_abstract_action_model(states, actions)
        opt_filename = os.path.join(config.experiment_dir, 'optimized.txt')
        logging.info("Minimized action model with {} actions saved in {}".format(len(actions), opt_filename))
        self.print_actions(actions, opt_filename, config.feature_namer)
        return states, actions, selected_data

    def print_actions(self, actions, filename, namer):
        with open(filename, 'w') as f:
            for i, action in enumerate(actions, 1):
                action_str = self.print_abstract_action(action, namer)
                print("\nAction {}:\n{}".format(i, action_str), file=f)

    def print_abstract_action(self, action, namer=lambda s: s):
        precs = ", ".join(sorted(self.print_precondition_atom(f, v, namer) for f, v in action.preconditions))
        effs = ", ".join(sorted(eff.print_named(namer) for eff in action.effects))
        return "\tPRE: {}\n\tEFFS: {}".format(precs, effs)

    def print_qnp_precondition_atom(self, feature, value, namer):
        assert value in (True, False)
        return "{} {}".format(namer(self.feature_names[feature]), int(value))

    def print_precondition_atom(self, feature, value, namer):
        assert value in (True, False)
        type_ = self.ftype(feature)
        fstr = namer(self.feature_names[feature])
        if type_ == bool:
            return fstr if value else "NOT {}".format(fstr)
        assert type_ == int
        return "{} > 0".format(fstr) if value else "{} = 0".format(fstr)

    def print_qnp_conjunction(self, conjunction, namer):
        sorted_ = sorted(self.print_qnp_precondition_atom(f, v, namer) for f, v in conjunction)
        return "{} {}".format(len(sorted_), " ".join(sorted_))

    def compute_qnp(self, states, actions, features, config, data):
        logging.info("Writing QNP abstraction to {}".format(config.qnp_abstraction_filename))
        namer = lambda x: config.feature_namer(x).replace(" ", "")  # let's make sure no spaces in feature names

        init_dnf, goal_dnf = self.compute_init_goals(data.sample, features)

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
        return self.minimize_dnf(dnf)

    def compute_init_goals(self, sample, features):
        fidxs = [f["idx"] for f in features]

        init_dnf = self.extract_dnf(fidxs, sorted(sample.roots))
        goal_dnf = self.extract_dnf(fidxs, sorted(sample.goals))

        return init_dnf, goal_dnf

    def minimize_dnf(self, dnf):
        return merge_precondition_sets(dnf)

    def find_inconsistency(self, features):
        abstract_states, abstraction = self.compute_state_abstractions(features)

        abstract2concrete = defaultdict(set)
        for s in self.optimal_and_nonoptimal:
            abstract2concrete[abstraction[s]].add(s)

        for s, sprime in self.iterate_over_optimal_transitions():
            abstract_s = abstraction[s]
            qs = self.retrieve_possibly_cached_qchanges(s, sprime)

            same_abstraction = [x for x in abstract2concrete[abstract_s] if x != s]

            for t in same_abstraction:
                if not self.sample.transitions[t]:  # i.e. state was not expanded
                    continue

                equal_qchanges_found = False
                all_qts = []
                for tprime in self.sample.transitions[t]:
                    qt = self.retrieve_possibly_cached_qchanges(t, tprime)
                    all_qts.append((tprime, qt))
                    if np.array_equal(qs, qt):
                        equal_qchanges_found = True
                        break
                if not equal_qchanges_found:
                    named_qs = list(zip(self.feature_names, qs))
                    logging.error("Unsoundness! No correspondence between qual. changes in transition {} "
                                  "and qual. changes of any transition starting at state {}".format((s, sprime), t))
                    logging.error("Q. changes for transition {}:".format((s, sprime)))
                    logging.error("\t{}".format(named_qs))
                    logging.error("Q. changes for all transitions starting at {}:".format(t))
                    for tprime, qt in all_qts:
                        named_qt = list(zip(self.feature_names, qt))
                        logging.error("\t{}: {}".format((t, tprime), named_qt))


def attempt_single_merge(action_precs):
    for p1, p2 in itertools.combinations(action_precs, 2):
        diff = p1.symmetric_difference(p2)
        diffl = list(diff)
        if len(diffl) == 2 and diffl[0][0] == diffl[1][0]:
            # The two conjunctions differ in that one has one literal L and the other its negation, the rest being equal
            assert diffl[0][1] != diffl[1][1]
            p_merged = p1.difference(diff)
            return p1, p2, p_merged  # Meaning p1 and p2 should be merged into p_merged
    return None


def merge_precondition_sets(action_precs):
    action_precs = action_precs.copy()
    while True:
        res = attempt_single_merge(action_precs)
        if res is None:
            break

        # else do the actual merge
        p1, p2, new = res
        action_precs.remove(p1)
        action_precs.remove(p2)
        action_precs.add(new)
    return action_precs


def optimize_abstract_action_model(states, actions):
    actions_grouped_by_effects = defaultdict(set)
    for a in actions:
        actions_grouped_by_effects[a.effects].add(a.preconditions)
    logging.info("Optimizing abstract action model with {} actions and {} effect groups".
                 format(len(actions), len(actions_grouped_by_effects)))

    merged = []
    for effs, action_precs in actions_grouped_by_effects.items():
        for prec in merge_precondition_sets(action_precs):
            merged.append(AbstractAction(prec, effs))

    return states, merged
