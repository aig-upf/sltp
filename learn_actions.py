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

from tarski.dl.features import relax_int_feature_diff

from errors import CriticalPipelineError
from tarski.dl import FeatureValueChange, MinDistanceFeature, NullaryAtomFeature, EmpiricalBinaryConcept, \
    ConceptCardinalityFeature
from util.console import print_header, print_lines
from util.command import count_file_lines, remove_duplicate_lines, read_file
from solvers import solve

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


def compute_qualitative_changes(transitions, all_features, model, relax_numeric_increase):
    # For each feature and each transition, precompute the qualitative change that the transition induces on the feature
    qchanges = {}
    for s0 in transitions:
        for s1 in transitions[s0]:
            for f in all_features:
                x0 = model.compute_feature_value(f, s0)
                x1 = model.compute_feature_value(f, s1)
                diff = f.diff(x0, x1)

                # If requested, we consider INCs and NILs as the same qualitative change.
                if relax_numeric_increase and isinstance(f, (MinDistanceFeature, )):
                    diff = relax_int_feature_diff(diff)
                qchanges[(s0, s1, f)] = diff
                # print("Denotation of feature \"{}\" along transition "
                #       "({},{}) is: {}".format(f, s0, s1, qchanges[(s0, s1, f)]))

    return qchanges


def compute_d1_distinguishing_features(states, transitions, features, model):
    ns, nf = len(states), len(features)
    logging.info("Computing sets of distinguishing features for each state pair over a total of "
                 "{} states and {} features ({:0.1f}K matrix entries)".format(ns, nf, nf*(ns*(ns-1))/(2*1000)))
    # self.undistinguishable_state_pairs = []
    d1_distinguishing_features = dict()

    # representative[s] will hold a pointer to the lowest-id state s' such that it is not possible to distinguish
    # s from s' with any of the full set of features
    for s1, s2 in itertools.combinations(states, 2):
        distinguishing = set()

        for f in features:
            x1 = model.compute_feature_value(f, s1)
            x2 = model.compute_feature_value(f, s2)
            if f.bool_value(x1) != f.bool_value(x2):  # f distinguishes s1 and s2
                distinguishing.add(f)

        d1_distinguishing_features[(s1, s2)] = distinguishing

    return states, transitions, d1_distinguishing_features

    # def transition_is_distinguishable(s1, t1, s2):
    #     for t2 in self.transitions[s2]:
    #         both_transitions_qualitatively_equal_over_all_features = True
    #         for f in self.features:
    #             change1 = self.qchanges[(s1, t1, f)]
    #             change2 = self.qchanges[(s2, t2, f)]
    #             if change1 != change2:
    #                 both_transitions_qualitatively_equal_over_all_features = False
    #                 break
    #         if not both_transitions_qualitatively_equal_over_all_features:
    #             # No feature in the global set of features is able to distinguish transition (s1,t1)
    #             # from some transition starting in s2
    #             can_ignore_s1 = False
    #
    #
    # for s1, s2 in representative.items():
    #     assert s1 > s2
    #     can_ignore_s1 = True
    #
    #     # Is there any transition starting in s1 such that for some feature in the global feature set
    #     # the transition can be qualitatively distinguished from all transitions starting in s2?
    #     # If not, then we can ignore s1 to all effects
    #
    #
    #
    #     for t1 in self.transitions[s1]:
    #         if transition_is_distinguishable(s1, t1, s2)
    #
    #         for t2 in self.transitions[s2]:
    #             both_transitions_qualitatively_equal_over_all_features = True
    #             for f in self.features:
    #                 change1 = self.qchanges[(s1, t1, f)]
    #                 change2 = self.qchanges[(s2, t2, f)]
    #                 if change1 != change2:
    #                     both_transitions_qualitatively_equal_over_all_features = False
    #                     break
    #             if not both_transitions_qualitatively_equal_over_all_features:
    #                 # No feature in the global set of features is able to distinguish transition (s1,t1)
    #                 # from some transition starting in s2
    #                 can_ignore_s1 = False


def generate_maxsat_problem(config, data):
    logging.info("Generating MAXSAT problem from {} concept-based features".format(len(data.features)))
    optimization = config.optimization if hasattr(config, "optimization") else OptimizationPolicy.NUM_FEATURES

    state_ids = data.state_ids
    features = data.features
    model = data.feature_model
    transitions = data.transitions
    goal_states = data.goal_states

    feat_matrix = np.load(config.feature_matrix_filename + ".npy")
    bin_feat_matrix = np.load(config.bin_feature_matrix_filename + ".npy")

    # Compute for each pair s and t of states which features distinguish s and t
    state_ids, transitions, d1_distinguishing_features = \
        compute_d1_distinguishing_features(state_ids, transitions, features, model)

    if not goal_states:
        raise CriticalPipelineError("No goal state identified in the sample, SAT theory will be trivially satisfiable")

    translator = ModelTranslator(features, state_ids, goal_states, transitions, model,
                                 config.cnf_filename, d1_distinguishing_features, optimization)

    translator.run(config.relax_numeric_increase)
    # translator.writer.print_variables(config.maxsat_variables_file)

    return dict(cnf_translator=translator)


def run_solver(config, data):
    # solution = solve(config.experiment_dir, config.cnf_filename, 'wpm3')
    # solution = solve(config.experiment_dir, config.cnf_filename, 'maxino')
    solution = solve(config.experiment_dir, config.cnf_filename, 'openwbo')
    if not solution.solved and solution.result == "UNSATISFIABLE":
        raise CriticalPipelineError("MAXSAT encoding is UNSATISFIABLE")
    else:
        logging.info("MAXSAT solution with cost {} found".format(solution.cost))

    return dict(cnf_translator=data.cnf_translator, cnf_solution=solution)


class AbstractAction(object):
    def __init__(self, name, preconditions, effects):
        self.name = name
        self.preconditions = frozenset(preconditions)
        self.effects = frozenset(effects)
        self.hash = hash((self.__class__, self.preconditions, self.effects))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.preconditions == other.preconditions and self.effects == other.effects)

    def __str__(self):
        return self.print_named()

    def print_named(self, namer=lambda s: s):
        precs = ", ".join(sorted(prec.print_named(namer) for prec in self.preconditions))
        effs = ", ".join(sorted(eff.print_named(namer) for eff in self.effects))
        return "\tPRE: {}\n\tEFFS: {}".format(precs, effs)


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

    def print_named(self, namer=lambda s: s):
        if self.change == FeatureValueChange.ADD:
            return namer(self.feature)
        if self.change == FeatureValueChange.DEL:
            return "NOT {}".format(namer(self.feature))
        if self.change == FeatureValueChange.INC:
            return "INC {}".format(namer(self.feature))
        if self.change == FeatureValueChange.DEC:
            return "DEC {}".format(namer(self.feature))
        if self.change == FeatureValueChange.INC_OR_NIL:
            return "INC* {}".format(namer(self.feature))
        raise RuntimeError("Unexpected effect type")


class AbstractPrecondition(Enum):
    ADD = 1
    DEL = 2
    GT0 = 3
    EQ0 = 4


def create_precondition_atom(feature, value):
    if value is True:
        return Atom(feature, AbstractPrecondition.ADD)
    if value is False:
        return Atom(feature, AbstractPrecondition.DEL)

    assert isinstance(value, int) and value >= 0

    if value > 0:
        return Atom(feature, AbstractPrecondition.GT0)

    return Atom(feature, AbstractPrecondition.EQ0)


class Atom(object):
    def __init__(self, feature, ptype):
        assert isinstance(ptype, AbstractPrecondition)
        self.feature = feature
        self.ptype = ptype
        self.hash = hash((self.__class__, self.feature, self.ptype))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.feature == other.feature and self.ptype == other.ptype)

    def __str__(self):
        return self.print_named()

    def print_named(self, namer=lambda s: s):
        if self.ptype is AbstractPrecondition.ADD:
            return namer(self.feature)
        if self.ptype is AbstractPrecondition.DEL:
            return "NOT {}".format(namer(self.feature))
        if self.ptype is AbstractPrecondition.GT0:
            return "{} > 0".format(namer(self.feature))
        if self.ptype is AbstractPrecondition.EQ0:
            return "{} = 0".format(namer(self.feature))
        assert False, "Unknown precondition type"


class ModelTranslator(object):
    def __init__(self, features, state_ids, goal_states, transitions, model, cnf_filename, d1_distinguishing_features,
                 optimization):
        self.state_ids = state_ids
        self.transitions = transitions
        self.model = model
        self.features = features
        self.goal_states = goal_states
        self.d1_distinguishing_features = d1_distinguishing_features

        self.writer = CNFWriter(cnf_filename)

        self.var_selected = None
        self.var_d1 = None
        self.var_d2 = None

        self.n_selected_clauses = 0
        self.n_d1_clauses = 0
        self.n_d2_clauses = 0
        self.n_bridge_clauses = 0
        self.n_goal_clauses = 0
        self.qchanges = None

        self.compute_feature_weight = self.setup_optimization_policy(optimization)

    def create_bridge_clauses(self, d1_literal, s, t):
        # If there are no transitions (t, t') in the sample set, then we do not post the bridge constraint.
        # if t not in self.transitions or not self.transitions[t]:
        if not self.transitions[t]:
            return

        for s_prime in self.transitions[s]:  # will be empty set if not initialized, which is ok
            forward_clause_literals = [d1_literal]
            for t_prime in self.transitions[t]:
                idx = compute_d2_index(s, s_prime, t, t_prime)
                forward_clause_literals.append(self.writer.literal(self.var_d2[idx], False))
            self.writer.clause(forward_clause_literals)
            self.n_bridge_clauses += 1

    def create_variables(self):
        print("Creating model variables".format())
        self.var_selected = {feat: self.writer.variable("selected({})".format(feat)) for feat in self.features}

        self.var_d1 = {(s1, s2): self.writer.variable("D1[{}, {}]".format(s1, s2)) for s1, s2 in
                       itertools.combinations(self.state_ids, 2)}

        self.var_d2 = dict()

        for s, t in itertools.combinations(self.transitions, 2):
            for s_prime, t_prime in itertools.product(self.transitions[s], self.transitions[t]):
                idx = compute_d2_index(s, s_prime, t, t_prime)
                varname = "D2[({},{}), ({},{})]".format(*idx)
                self.var_d2[idx] = self.writer.variable(varname)

        self.writer.close()
        print("A total of {} variables were created".format(len(self.writer.variables)))

    def run(self, relax_numeric_increase):
        self.qchanges = qchanges = compute_qualitative_changes(self.transitions, self.features, self.model,
                                                               relax_numeric_increase)

        self.create_variables()

        print("Generating D1 + bridge constraints from {} D1 variables".format(len(self.var_d1)), flush=True)
        for s1, s2 in itertools.combinations(self.state_ids, 2):
            d1_variable = self.var_d1[(s1, s2)]
            d1_distinguishing_features = self.d1_distinguishing_features[(s1, s2)]

            # Post the constraint: D1(si, sj) <=> OR active(f), where the OR ranges over all
            # those features f that tell apart si from sj
            d1_lit = self.writer.literal(d1_variable, True)
            forward_clause_literals = [self.writer.literal(d1_variable, False)]

            for f in d1_distinguishing_features:
                forward_clause_literals.append(self.writer.literal(self.var_selected[f], True))
                self.writer.clause([d1_lit, self.writer.literal(self.var_selected[f], False)])
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
                self.create_bridge_clauses(d1_lit, s1, s2)
                self.create_bridge_clauses(d1_lit, s2, s1)
                pass

        print("Generating D2 constraints from {} D2 variables".format(len(self.var_d2)), flush=True)
        # self.d2_distinguished = set()
        for (s0, s1, t0, t1), d2_var in self.var_d2.items():
            d2_distinguishing_features = []  # all features that d2-distinguish the current pair of transitions
            for f in self.features:
                if qchanges[(s0, s1, f)] != qchanges[(t0, t1, f)]:
                    d2_distinguishing_features.append(f)

            # if len(d2_distinguishing_features) != 0:
            #     self.d2_distinguished.add((s0, s1, t0, t1))

            # D2(s0,s1,t0,t2) iff OR_f selected(f), where f ranges over features that d2-distinguish the transition
            # but do _not_ d1-distinguish the two states at the origin of each transition.
            d2_lit = self.writer.literal(d2_var, True)
            forward_clause_literals = [self.writer.literal(d2_var, False)]
            # for f in (x for x in d2_distinguishing_features):
            for f in (x for x in d2_distinguishing_features if x not in self.d1_distinguishing_features[(s0, t0)]):
                forward_clause_literals.append(self.writer.literal(self.var_selected[f], True))
                self.writer.clause([d2_lit, self.writer.literal(self.var_selected[f], False)])
                self.n_d2_clauses += 1

            self.writer.clause(forward_clause_literals)
            self.n_d2_clauses += 1

        # Add the weighted clauses to minimize the number of selected features
        for feature, feat_var in self.var_selected.items():
            d = self.compute_feature_weight(feature)
            self.writer.clause([self.writer.literal(feat_var, False)], weight=d)
            self.n_selected_clauses += 1

        self.report_stats()

        # self.debug_tests()

        self.writer.save()

        return self.writer.mapping

    def decode_solution(self, assignment, namer):
        varmapping = self.writer.mapping
        true_variables = set(varmapping[idx] for idx, value in assignment.items() if value is True)
        feature_mapping = {variable: feature for feature, variable in self.var_selected.items()}
        assert len(feature_mapping) == len(self.var_selected)
        selected_features = [feature_mapping[v] for v in true_variables if v in feature_mapping]

        if not selected_features:
            raise CriticalPipelineError("Zero-cost maxsat solution - "
                                        "no action model possible, the encoding has likely some error")

        logging.info("Selected features: ")
        print('\n'.join("F{}. {} [{}]".format(i, namer(f), f.complexity()) for i, f in enumerate(selected_features, 1)))

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

                if (abstract_s, abstract_sprime) in already_computed:
                    continue

                abstract_effects = []
                for f in selected_features:
                    concrete_qchange = self.qchanges[(s, sprime, f)]
                    if concrete_qchange != FeatureValueChange.NIL:
                        abstract_effects.append(ActionEffect(f, concrete_qchange))

                preconditions = [create_precondition_atom(f, val) for f, val in zip(selected_features, abstract_s)]
                # preconditions = [Atom(f, val) for f, val in zip(selected_features, abstract_s)]
                abstract_actions.add(AbstractAction("anonymous", preconditions,
                                                    abstract_effects))
                if len(abstract_effects) == 0:
                    raise RuntimeError("Unsound state model abstraction!")

                already_computed.add((abstract_s, abstract_sprime))

        logging.info("Abstract state space: {} states and {} actions".
                     format(len(abstract_states), len(abstract_actions)))
        return abstract_states, abstract_actions

    def report_stats(self):

        d1_distinguishing = self.d1_distinguishing_features.values()
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
        return feature.complexity()

    def setup_optimization_policy(self, optimization):
        if optimization == OptimizationPolicy.NUM_FEATURES:
            return self.opt_policy_num_features

        if optimization == OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY:
            return self.opt_policy_total_feature_depth

        raise RuntimeError("Unknown optimization policy")

    def compute_abstract_state(self, features, state_id):
        state_values = []
        for f in features:
            x = self.model.compute_feature_value(f, state_id)
            state_values.append(x)

        # if f.bool_value(x1) != f.bool_value(x2):  # f distinguishes s1 and s2

        return tuple(state_values)


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


def compute_action_model(config, data):
    assert data.cnf_solution.solved

    states, actions = data.cnf_translator.decode_solution(data.cnf_solution.assignment, config.feature_namer)
    with open(os.path.join(config.experiment_dir, 'actions.txt'), 'w') as f:
        f.write("\n\n".join("Action {}:\n{}".format(i, action.print_named(config.feature_namer))
                            for i, action in enumerate(actions, 1)))
    return dict(abstract_states=states, abstract_actions=actions)
