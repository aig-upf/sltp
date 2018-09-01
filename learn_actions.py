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
from signal import signal, SIGPIPE, SIG_DFL
import time

import numpy as np

from errors import CriticalPipelineError
from tarski.dl import FeatureValueChange, MinDistanceFeature
from util.console import print_header, print_lines
from util.command import count_file_lines, remove_duplicate_lines, read_file
from solvers import solve
from util.performance import print_memory_usage
from util.serialization import serialize_to_string

signal(SIGPIPE, SIG_DFL)

BASEDIR = os.path.dirname(os.path.realpath(__file__))


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
    optimization = config.optimization if hasattr(config, "optimization") else OptimizationPolicy.NUM_FEATURES

    state_ids = data.state_ids
    optimal_state_ids = data.optimal_states
    transitions = data.transitions
    goal_states = data.goal_states

    feature_complexity = np.load(config.feature_complexity_filename + ".npy")
    feature_names = np.load(config.feature_names_filename + ".npy")
    feat_matrix = np.load(config.feature_matrix_filename + ".npy")
    bin_feat_matrix = np.load(config.bin_feature_matrix_filename + ".npy")
    feature_types = feat_matrix.max(0) <= 1  # i.e. feature_types[i] is True iff i is an (empirically) bool feature

    logging.info("Generating MAXSAT problem from {} states and {} features"
                 .format(feat_matrix.shape[0], feat_matrix.shape[1]))

    if not goal_states:
        raise CriticalPipelineError("No goal state identified in the sample, SAT theory will be trivially satisfiable")

    translator = ModelTranslator(feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types,
                                 state_ids, goal_states, transitions, optimal_state_ids,
                                 config.cnf_filename, optimization,
                                 config.relax_numeric_increase, config.complete_only_wrt_optimal)

    translator.run(data.enforced_feature_idxs)
    # translator.writer.print_variables(config.maxsat_variables_file)

    return dict(cnf_translator=translator)


def run_solver(config, data, rng):
    # solution = solve(config.experiment_dir, config.cnf_filename, 'wpm3')
    # solution = solve(config.experiment_dir, config.cnf_filename, 'maxino')
    solution = solve(config.experiment_dir, config.cnf_filename, 'openwbo')
    if not solution.solved and solution.result == "UNSATISFIABLE":
        raise CriticalPipelineError("MAXSAT encoding is UNSATISFIABLE")
    else:
        logging.info("MAXSAT solution with cost {} found".format(solution.cost))

    return dict(cnf_solution=solution)


class AbstractAction(object):
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
        return "AbstractAction[]"

    __repr__ = __str__


class ActionEffect(object):
    def __init__(self, feature, change):
        self.feature = feature
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
        if self.change == FeatureValueChange.ADD:
            return namer(self.feature)
        if self.change == FeatureValueChange.DEL:
            return "NOT {}".format(namer(self.feature))
        if self.change == FeatureValueChange.ADD_OR_NIL:
            return "ADD* {}".format(namer(self.feature))
        if self.change == FeatureValueChange.INC:
            return "INC {}".format(namer(self.feature))
        if self.change == FeatureValueChange.DEC:
            return "DEC {}".format(namer(self.feature))
        if self.change == FeatureValueChange.INC_OR_NIL:
            return "INC* {}".format(namer(self.feature))
        raise RuntimeError("Unexpected effect type")

    def print_qnp_named(self, namer=lambda s: s):
        name = namer(self.feature)
        if self.change in (FeatureValueChange.ADD, FeatureValueChange.INC):
            return "{} 1".format(name)

        if self.change in (FeatureValueChange.DEL, FeatureValueChange.DEC):
            return "{} 0".format(name)

        if self.change in (FeatureValueChange.INC_OR_NIL, FeatureValueChange.ADD_OR_NIL):
            assert False, "Relaxed INC semantics not supported for QNP"
        raise RuntimeError("Unexpected effect type")


class ModelTranslator(object):
    def __init__(self, feat_matrix, bin_feat_matrix, feature_complexity, feature_names, feature_types, state_ids,
                 goal_states, transitions, optimal_state_ids, cnf_filename,
                 optimization, relax_numeric_increase, complete_only_wrt_optimal):

        self.feat_matrix = feat_matrix
        self.bin_feat_matrix = bin_feat_matrix
        self.feature_complexity = feature_complexity
        self.feature_names = feature_names
        self.feature_types = feature_types
        self.state_ids = state_ids
        self.optimal_state_ids = optimal_state_ids
        self.transitions = transitions
        self.goal_states = goal_states
        self.complete_only_wrt_optimal = complete_only_wrt_optimal

        # We'll need this later if using relaxed inc semantics:
        self.all_twos = 2*np.ones(self.feat_matrix.shape[1], dtype=np.int8)

        # Compute for each pair s and t of states which features distinguish s and t
        if complete_only_wrt_optimal:
            self.optimal_states = optimal_state_ids
            self.optimal_and_nonoptimal = set(state_ids).union(optimal_state_ids)
        else:
            assert False
            # self.completeness_and_correctness_states = set(state_ids)
            # self.correctness_only_states = set()

        self.np_d1_distinguishing_features = self.compute_d1_distinguishing_features(bin_feat_matrix)

        self.writer = CNFWriter(cnf_filename)

        self.np_var_selected = None
        self.var_d1 = None
        self.var_d2 = None

        self.n_selected_clauses = 0
        self.n_d1_clauses = 0
        self.n_d2_clauses = 0
        self.n_bridge_clauses = 0
        self.n_goal_clauses = 0
        self.np_qchanges = dict()
        self.compute_qchanges = self.relaxed_qchange if relax_numeric_increase else self.standard_qchange
        self.compute_feature_weight = self.setup_optimization_policy(optimization)

    def iterate_over_d1_pairs(self):
        for (x, y) in itertools.product(self.optimal_states, self.optimal_and_nonoptimal):
            if y in self.optimal_states and y <= x:
                continue
            yield (x, y)

    def compute_d1_distinguishing_features(self, bin_feat_matrix):
        """ """
        # assert len(set(self.optimal_states) & set(self.optimal_and_nonoptimal)) == 0
        ns1, ns2, nf = len(self.optimal_states), len(self.optimal_and_nonoptimal), bin_feat_matrix.shape[1]

        # How many feature denotations we'll have to compute
        nentries = nf * ((ns1 * (ns1-1)) / 2 + ns1 * ns2)

        logging.info("Computing sets of D1-distinguishing features for {} x {} state pairs "
                     "and {} features ({:0.1f}M matrix entries)".format(ns1, ns2, nf, nentries / 1000000))
        np_d1_distinguishing_features = dict()

        for s1, s2 in self.iterate_over_d1_pairs():
            xor = np.logical_xor(bin_feat_matrix[s1], bin_feat_matrix[s2])
            # We extract the tuple and turn it into a set:
            np_d1_distinguishing_features[(s1, s2)] = set(np.nonzero(xor)[0].flat)

        return np_d1_distinguishing_features

    def create_bridge_clauses(self, d1_literal, s, t):
        # If there are no transitions (t, t') in the sample set, then we do not post the bridge constraint.
        # if t not in self.transitions or not self.transitions[t]:
        # TODO WE MIGHT NEED A BETTER WAY OF IDENTIFYING NON-EXPANDED STATES AND DISTINGUISHING THEM FROM
        # TODO STATES WHICH HAVE NO SUCCESSOR
        if not self.transitions[s] or not self.transitions[t]:
            return

        assert s != t
        assert s in self.optimal_states

        for s_prime in self.transitions[s]:
            if s_prime < s or s_prime not in self.optimal_states:  # We want only optimal transitions
                continue

            forward_clause_literals = [d1_literal]
            for t_prime in self.transitions[t]:
                idx = compute_d2_index(s, s_prime, t, t_prime)
                forward_clause_literals.append(self.writer.literal(self.var_d2[idx], False))
            self.writer.clause(forward_clause_literals)
            self.n_bridge_clauses += 1

    def retrieve_possibly_cached_qchanges(self, s0, s1):
        try:
            return self.np_qchanges[(s0, s1)]
        except KeyError:
            self.np_qchanges[(s0, s1)] = val = self.compute_qchanges(s0, s1)
            return val

    def iterate_over_optimal_transitions(self):
        """ """
        for s in self.optimal_states:
            for sprime in self.transitions[s]:
                if sprime < s or sprime not in self.optimal_states:
                    continue
                yield (s, sprime)

    def iterate_over_all_transitions(self):
        for s in self.transitions:
            for s_prime in self.transitions[s]:
                yield (s, s_prime)

    def iterate_over_d2_transitions(self):
        for s, s_prime in self.iterate_over_optimal_transitions():
            for t, t_prime in self.iterate_over_all_transitions():
                if s == t or (self.is_optimal_transition(t, t_prime) and t < s):
                    continue
                yield (s, s_prime, t, t_prime)

    def is_optimal_transition(self, s, sprime):
        return s in self.optimal_states and sprime in self.optimal_states and s < sprime

    def create_variables(self):
        # self.var_selected = {feat: self.writer.variable("selected({})".format(feat)) for feat in self.features}
        self.np_var_selected = \
            [self.writer.variable("selected(F{})".format(feat)) for feat in range(0, self.feat_matrix.shape[1])]

        self.var_d1 = {(s1, s2): self.writer.variable("D1[{}, {}]".format(s1, s2)) for s1, s2 in
                       self.iterate_over_d1_pairs()}

        self.var_d2 = dict()

        for s, s_prime, t, t_prime in self.iterate_over_d2_transitions():
            # for s, t in itertools.combinations(self.transitions, 2):
            #     for s_prime, t_prime in itertools.product(self.transitions[s], self.transitions[t]):
            idx = compute_d2_index(s, s_prime, t, t_prime)
            varname = "D2[({},{}), ({},{})]".format(*idx)
            self.var_d2[idx] = self.writer.variable(varname)

        self.writer.close()
        logging.info("A total of {} MAXSAT variables were created".format(len(self.writer.variables)))

    def compute_d1_idx(self, s, t):
        if s in self.optimal_states and t in self.optimal_states:
            if s > t:
                return t, s
        return s, t

    def run(self, enforced_feature_idxs):

        self.create_variables()

        print("Generating D1 + bridge constraints from {} D1 variables".format(len(self.var_d1)), flush=True)

        # for s1, s2 in itertools.combinations(self.state_ids, 2):
        for s1, s2 in self.iterate_over_d1_pairs():
            d1_variable = self.var_d1[(s1, s2)]
            # d1_distinguishing_features = self.d1_distinguishing_features[(s1, s2)]
            d1_distinguishing_features = self.np_d1_distinguishing_features[(s1, s2)]

            # Post the constraint: D1(si, sj) <=> OR active(f), where the OR ranges over all
            # those features f that tell apart si from sj
            d1_lit = self.writer.literal(d1_variable, True)
            forward_clause_literals = [self.writer.literal(d1_variable, False)]

            for f in d1_distinguishing_features:
                forward_clause_literals.append(self.writer.literal(self.np_var_selected[f], True))
                self.writer.clause([d1_lit, self.writer.literal(self.np_var_selected[f], False)])
                self.n_d1_clauses += 1

            self.writer.clause(forward_clause_literals)
            self.n_d1_clauses += 1

            # Force D1(s1, s2) to be true if exactly one of the two states is a goal state
            if sum(1 for x in (s1, s2) if x in self.goal_states) == 1:
                self.writer.clause([d1_lit])
                self.n_goal_clauses += 1

                if len(d1_distinguishing_features) == 0:
                    print("WARNING: No feature in the pool able to distinguish state pair {}, but one of the states "
                          "is a goal and the other is not. MAXSAT encoding will be UNSAT".format((s1, s2)), flush=True)

            # Else (i.e. D1(s1, s2) _might_ be false, create the bridge clauses between values of D1 and D2
            else:
                assert s1 in self.optimal_states
                self.create_bridge_clauses(d1_lit, s1, s2)
                if s2 in self.optimal_states:
                    self.create_bridge_clauses(d1_lit, s2, s1)

        print("Generating D2 constraints from {} D2 variables".format(len(self.var_d2)), flush=True)
        # self.d2_distinguished = set()
        # for (s0, s1, t0, t1), d2_var in self.var_d2.items():
        for s0, s1, t0, t1 in self.iterate_over_d2_transitions():
            d2_var = self.var_d2[compute_d2_index(s0, s1, t0, t1)]

            qchanges_s0s1 = self.retrieve_possibly_cached_qchanges(s0, s1)
            qchanges_t0t1 = self.retrieve_possibly_cached_qchanges(t0, t1)
            equal_idxs = np.equal(qchanges_s0s1, qchanges_t0t1)
            np_d2_distinguishing_features = np.where(equal_idxs==False)[0]

            np_d1_dist = self.np_d1_distinguishing_features[self.compute_d1_idx(s0, t0)]
            # if s0 == 440:
            #     print(np_d1_dist)
            # D2(s0,s1,t0,t2) iff OR_f selected(f), where f ranges over features that d2-distinguish the transition
            # but do _not_ d1-distinguish the two states at the origin of each transition.
            d2_lit = self.writer.literal(d2_var, True)
            forward_clause_literals = [self.writer.literal(d2_var, False)]
            for f in np_d2_distinguishing_features.flat:
                if f not in np_d1_dist:
                    forward_clause_literals.append(self.writer.literal(self.np_var_selected[f], True))
                    self.writer.clause([d2_lit, self.writer.literal(self.np_var_selected[f], False)])
                    self.n_d2_clauses += 1

            self.writer.clause(forward_clause_literals)
            self.n_d2_clauses += 1

        # Add the weighted clauses to minimize the number of selected features
        for feature, var in enumerate(self.np_var_selected):
            d = self.compute_feature_weight(feature)
            self.writer.clause([self.writer.literal(var, False)], weight=d)
            self.n_selected_clauses += 1

        # Make sure those features we want to enforce in the solution, if any, are posted as necessarily selected
        for enforced in enforced_feature_idxs:
            self.writer.clause([self.writer.literal(self.np_var_selected[enforced], True)])

        # from pympler.asizeof import asizeof
        # print(f"asizeof(self.writer): {asizeof(self.writer)/(1024*1024)}MB")
        # print(f"asizeof(self): {asizeof(self)/(1024*1024)}MB")
        # print(f"asizeof(self.np_d1_distinguishing_features): {asizeof(self.np_d1_distinguishing_features)/(1024*1024)}MB")
        # print(f"np_d2_distinguishing_features.nbytes: {np_d2_distinguishing_features.nbytes/(1024*1024)}MB")
        print_memory_usage()

        self.report_stats()

        # self.debug_tests()
        self.np_d1_distinguishing_features = None  # We won't need this anymore
        self.writer.save()

        return self.writer.mapping

    def standard_qchange(self, s0, s1):
        return np.sign(self.feat_matrix[s1] - self.feat_matrix[s0])

    def relaxed_qchange(self, s0, s1):
        std = self.standard_qchange(s0, s1)
        # For each feature, return the standard qchange value if the feature is boolean, or the qchange is -1 (DEC),
        # but otherwise turn 0's (NO-CHANGE) and 1's (INCREASES) into value 2, denoted to mean "relaxed increase" (INC*)
        relaxed = np.where(np.logical_or(std == -1, self.feature_types == True), std, self.all_twos)
        return relaxed

    def ftype(self, f):
        return bool if self.feature_types[f] else int

    def generate_eff_change(self, feat, qchange):

        if self.ftype(feat) == bool:
            assert qchange in (-1, 1)
            return {-1: FeatureValueChange.DEL, 1: FeatureValueChange.ADD, 2: FeatureValueChange.ADD_OR_NIL}[qchange]

        # else type must be int
        assert qchange in (-1, 1, 2)
        return {-1: FeatureValueChange.DEC, 1: FeatureValueChange.INC, 2: FeatureValueChange.INC_OR_NIL}[qchange]

    def decode_solution(self, assignment, features, namer):
        varmapping = self.writer.mapping
        true_variables = set(varmapping[idx] for idx, value in assignment.items() if value is True)
        # feature_mapping = {variable: feature for feature, variable in self.var_selected.items()}
        np_feature_mapping = {variable: feature for feature, variable in enumerate(self.np_var_selected)}
        assert len(np_feature_mapping) == len(self.np_var_selected)
        selected_features = sorted(np_feature_mapping[v] for v in true_variables if v in np_feature_mapping)

        if not selected_features:
            raise CriticalPipelineError("Zero-cost maxsat solution - "
                                        "no action model possible, the encoding has likely some error")

        logging.info("Features (total complexity: {}): ".format(sum(self.feature_complexity[f] for f in selected_features)))
        print('\t' + '\n\t'.join("{}. {} [k={}, id={}]".format(i, namer(self.feature_names[f]), self.feature_complexity[f], f)
                        for i, f in enumerate(selected_features, 1)))

        print("Only for machine eyes: ")
        print(serialize_to_string([features[i] for i in selected_features]))
        selected_data = [dict(idx=f,
                              name=namer(self.feature_names[f]),
                              type=int(not self.feature_types[f]))
                         for f in selected_features]

        # selected_features = [f for f in selected_features if str(f) == "bool[And(clear,{a})]"]

        abstract_states = set()
        state_abstraction = dict()
        for state in self.state_ids:
            sprime = self.compute_abstract_state(selected_features, state)
            abstract_states.add(sprime)
            state_abstraction[state] = sprime

        abstract_actions = set()
        already_computed = set()
        for s, children in self.transitions.items():
            for sprime in children:
                abstract_s = state_abstraction[s]
                abstract_sprime = state_abstraction[sprime]

                qchanges = self.retrieve_possibly_cached_qchanges(s, sprime)
                selected_qchanges = qchanges[selected_features]
                abstract_effects = [ActionEffect(self.feature_names[f], self.generate_eff_change(f, c))
                                    for f, c in zip(selected_features, selected_qchanges) if c != 0]

                precondition_bitmap = frozenset(zip(selected_features, abstract_s))
                abstract_actions.add(AbstractAction(precondition_bitmap, abstract_effects))
                if len(abstract_effects) == 0:
                    msg = "Abstract no-op necessary [concrete: ({}, {}), abstract: ({}, {})]"\
                        .format(s, sprime, abstract_s, abstract_sprime)
                    logging.warning(msg)
                    # raise RuntimeError(msg)

                already_computed.add((abstract_s, abstract_sprime))

        logging.info("Abstract state space: {} states and {} actions".
                     format(len(abstract_states), len(abstract_actions)))
        return abstract_states, abstract_actions, selected_data

    def report_stats(self):

        d1_distinguishing = self.np_d1_distinguishing_features.values()
        avg_num_d1_dist_features = sum(len(feats) for feats in d1_distinguishing) / len(d1_distinguishing)
        n_undistinguishable_state_pairs = sum(1 for x in d1_distinguishing if len(x) == 0)
        print_header("Max-sat encoding stats", 1)
        print_lines("Number of D1-undistinguishable state pairs: {}".format(n_undistinguishable_state_pairs), 1)
        print_lines("Avg. # of D1-distinguishing features: {:.2f}".format(avg_num_d1_dist_features), 1)
        print_lines("Clauses (possibly with repetitions):".format(), 1)
        print_lines("Selected: {}".format(self.n_selected_clauses), 2)
        print_lines("D1: {}".format(self.n_d1_clauses), 2)
        print_lines("D2: {}".format(self.n_d2_clauses), 2)
        print_lines("Bridge: {}".format(self.n_bridge_clauses), 2)
        print_lines("Goal: {}".format(self.n_goal_clauses), 2)
        print_lines("TOTAL: {}".format(self.n_selected_clauses + self.n_d1_clauses + self.n_d2_clauses +
                                       self.n_bridge_clauses + self.n_goal_clauses), 2)

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
        # state_values = []
        # for f in features:
        #     x = self.model.compute_feature_value(f, state_id)
        #     state_values.append(x)

        bin_values = self.bin_feat_matrix[state_id, features]
        return tuple(map(bool, bin_values))

    def compute_action_model(self, assignment, features, config):
        states, actions, selected_features = self.decode_solution(assignment, features, config.feature_namer)
        self.print_actions(actions, os.path.join(config.experiment_dir, 'actions.txt'), config.feature_namer)
        states, actions = optimize_abstract_action_model(states, actions)
        self.print_actions(actions, os.path.join(config.experiment_dir, 'actions-optimized.txt'), config.feature_namer)
        return states, actions, selected_features

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

    def compute_qnp(self, states, actions, features, config, data):
        logging.info("Writing QNP abstraction to {}".format(config.qnp_abstraction_filename))
        namer = lambda x: config.feature_namer(x).replace(" ", "")  # let's make sure no spaces in feature names

        init, goal = self.compute_init_goals(data, features)

        init = next(iter(init))  # TODO ATM we just pick the first initial state, without minimization
        goal = next(iter(goal))  # TODO ATM we just pick the first goal state, without minimization

        with open(config.qnp_abstraction_filename, "w") as f:
            print(config.domain_dir, file=f)

            name_types = " ".join("{} {}".format(namer(f["name"]), f["type"]) for f in features)
            feat_line = "{} {}".format(len(features), name_types)
            print(feat_line, file=f)

            init_conds = [self.print_qnp_precondition_atom(f, v, namer) for f, v in init]
            goal_conds = [self.print_qnp_precondition_atom(f, v, namer) for f, v in goal]
            print("{} {}".format(len(init_conds), " ".join(init_conds)), file=f)
            print("{} {}".format(len(goal_conds), " ".join(goal_conds)), file=f)

            # Actions
            print("{}".format(len(actions)), file=f)
            for i, action in enumerate(actions, 1):
                precs = sorted(self.print_qnp_precondition_atom(f, v, namer) for f, v in action.preconditions)
                effs = sorted(eff.print_qnp_named(namer) for eff in action.effects)
                print("action_{}".format(i), file=f)
                print("{} {}".format(len(precs), " ".join(precs)), file=f)
                print("{} {}".format(len(effs), " ".join(effs)), file=f)

    def extract_dnf(self, feature_indexes, states):
        """ """
        dnf = self.bin_feat_matrix[np.ix_(states, feature_indexes)]  # see https://stackoverflow.com/q/22927181
        dnf = set(frozenset((i, val) for i, val in zip(feature_indexes, f)) for f in dnf)
        return self.minimize_dnf(dnf)

    def compute_init_goals(self, data, features):
        fidxs = [f["idx"] for f in features]

        init_dnf = self.extract_dnf(fidxs, sorted(data.root_states))
        goal_dnf = self.extract_dnf(fidxs, sorted(data.goal_states))

        return init_dnf, goal_dnf

    def minimize_dnf(self, dnf):
        return merge_precondition_sets(dnf)


class Variable(object):
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name

    __repr__ = __str__


class Literal(object):
    def __init__(self, variable, polarity=True):
        assert isinstance(variable, Variable)
        self.variable = variable
        self.polarity = polarity

    def __neg__(self):
        return Literal(self.variable, not self.polarity)

    def __str__(self):
        return "{}{}".format(("" if self.polarity else "-"), self.variable)

    def to_cnf(self, variable_index):
        return "{}{}".format(("" if self.polarity else "-"), variable_index[self.variable])

    def __eq__(self, other):
        return self.variable == other.variable and self.polarity == other.polarity

    def __hash__(self):
        return hash((self.variable, self.polarity))

    __repr__ = __str__


class Clause(object):
    def __init__(self, literals, weight=None):
        # assert all(isinstance(l, Literal) for l in literals)
        # self.literals = tuple(literals)
        self.literals = literals
        # assert len(set(literals)) == len(self.literals)  # Make sure all literals are unique
        self.weight = weight

    def __str__(self):
        return "{{{}}} [{}]".format(','.join(str(l) for l in self.literals), self.weight)

    __repr__ = __str__

    # def cnf_line(self, top, variable_index):
    #     # w <literals> 0
    #     w = top if self.weight is math.inf else self.weight
    #     literals = " ".join(l.to_cnf(variable_index) for l in self.literals)
    #     return "{} {} 0".format(w, literals)

    # def cnf_line(self, top):
    #     # w <literals> 0
    #     w = top if self.weight is math.inf else self.weight
    #     literals = " ".join(str(l) for l in self.literals)
    #     return "{} {} 0".format(w, literals)

    # def __eq__(self, other):
    #     return self.literals == other.literals
    #
    # def __hash__(self):
    #     return hash(self.literals)


def print_clause(literals, weight):
    # w <literals> 0
    w = "#TOP#" if weight is None else weight
    literals = " ".join(str(l) for l in literals)
    return "{} {} 0".format(w, literals)


class CNFWriter(object):
    def __init__(self, filename):
        self.filename = filename
        self.variables = dict()
        self.clauses = set()
        self.clause_batch = []
        self.num_clauses = 0
        self.accumulated_weight = 0
        self.mapping = dict()
        self.closed = False
        self.variable_index = None

        self.buffer_filename = self.filename + ".tmp"
        self.buffer = open(self.buffer_filename, "w")

    def variable(self, name):
        assert not self.closed
        return self.variables.setdefault(name, Variable(name))

    def literal(self, variable, polarity):
        assert self.closed
        i = self.variable_index[variable]
        return i if polarity else -1 * i

    def clause(self, literals, weight=None):
        assert self.closed
        # self.clauses.add(Clause(literals=literals, weight=weight)) # Keeping the set in memory is expensive!
        self.num_clauses += 1
        self.accumulated_weight += weight if weight is not None else 0

        self.clause_batch.append((literals, weight))
        if len(self.clause_batch) == 1000:
            self.flush_clauses()

    def flush_clauses(self):
        clauses_str = '\n'.join(print_clause(literals, weight) for literals, weight in self.clause_batch)
        print(clauses_str, file=self.buffer)
        del self.clause_batch
        self.clause_batch = []

    def save(self):
        assert self.closed
        numvars = len(self.variables)
        numclauses = self.num_clauses

        self.flush_clauses()
        self.buffer.close()
        self.buffer = None
        self._save(self.filename, numvars, numclauses)

        # debug = True
        # debug = False
        # if debug:
        #     dfilename = "{}.txt".format(self.filename)
        #     debug_clause_printer = lambda c: str(c)
        #     self._save(dfilename, numvars, numclauses, top, debug_clause_printer)

    def _save(self, filename, numvars, numclauses):
        num_written_clauses = count_file_lines(self.buffer_filename)
        assert numclauses == num_written_clauses
        remove_duplicate_lines(self.buffer_filename)
        num_unique_clauses = count_file_lines(self.buffer_filename)
        # num_unique_clauses_in_mem = len(self.clauses)
        # assert num_unique_clauses == num_unique_clauses_in_mem  # Keeping the set in memory is expensive!
        top = str(self.accumulated_weight + 1)

        print("Writing max-sat encoding to file \"{}\"".format(self.filename))
        print("Max-sat problem has {} variables and {} unique clauses (with repetitions: {}). Top weight is {}".format(
            numvars, num_unique_clauses, numclauses, top))

        with open(filename, "w") as output:
            print("c WCNF model generated on {}".format(time.strftime("%Y%m%d %H:%M:%S", time.localtime())),
                  file=output)
            # p wcnf nbvar nbclauses top
            print("p wcnf {} {} {}".format(numvars, num_unique_clauses, top), file=output)
            for line in read_file(self.buffer_filename):
                print(line.replace("#TOP#", top), file=output)
                # for clause in self.clauses:
                #     print(clause_printer(clause), file=file)

    def close(self):  # Once closed, the writer won't admit more variable declarations
        assert not self.closed
        self.closed = True
        self.variable_index = {var: i for i, var in enumerate(self.variables.values(), start=1)}
        # Save the variable mapping to parse the solution later
        self.mapping = {i: name for name, i in self.variable_index.items()}

    def print_variables(self, filename):
        assert self.closed
        variables = sorted(self.variable_index.items(), key=lambda x: x[1])  # Sort by var index
        with open(filename, 'w') as f:
            print("\n".join("{}: {}".format(i, v) for v, i in variables), file=f)


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
    logging.info("Optimizing abstract action model with {} abstract actions".format(len(actions)))
    actions_grouped_by_effects = defaultdict(set)
    for a in actions:
        actions_grouped_by_effects[a.effects].add(a.preconditions)
    logging.info("There are {} effect groups".format(len(actions_grouped_by_effects)))

    merged = []
    for effs, action_precs in actions_grouped_by_effects.items():
        for prec in merge_precondition_sets(action_precs):
            merged.append(AbstractAction(prec, effs))

    logging.info("Optimized abstract action model has {} actions".format(len(merged)))
    return states, merged


def compute_action_model(config, data, rng):
    assert data.cnf_solution.solved
    states, actions, features = data.cnf_translator.compute_action_model(data.cnf_solution.assignment, data.features, config)
    data.cnf_translator.compute_qnp(states, actions, features, config, data)
    # return dict(abstract_states=states, abstract_actions=actions)
    return dict()
