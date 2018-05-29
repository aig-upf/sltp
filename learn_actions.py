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
import math
import os
from collections import defaultdict
from signal import signal, SIGPIPE, SIG_DFL
import time

from extensions import ExtensionCache
from util.console import print_header, print_lines
from util.command import count_file_lines, remove_duplicate_lines, read_file
from solvers import solve

signal(SIGPIPE, SIG_DFL)

BASEDIR = os.path.dirname(os.path.realpath(__file__))


def compute_d2_index(s0, s1, t0, t1):
    assert s0 != t0
    if s0 < t0:
        return s0, s1, t0, t1
    else:
        return t0, t1, s0, s1


def compute_qualitative_changes(transitions, all_features, model, substitution):
    # For each feature and each transition, precompute the qualitative change that the transition induces on the feature
    qchanges = {}
    for s0 in transitions:
        for s1 in transitions[s0]:
            for f in all_features:
                x0 = model.compute_feature_value(f, s0, substitution)
                x1 = model.compute_feature_value(f, s1, substitution)
                qchanges[(s0, s1, f)] = f.diff(x0, x1)
                # print("Denotation of feature \"{}\" along transition ({},{}) is: {}".format(f, s0, s1, qchanges[(s0, s1, f)]))

    return qchanges


def run(config, data):
    translator = ModelTranslator(data.features, data.states, data.goal_states, data.transitions, data.extensions,
                                 config.cnf_filename)
    translator.log_features(config.feature_filename)

    # DEBUGGING: FORCE A CERTAIN SET OF FEATURES:
    selected = None
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 1232]]
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 2390]]
    # selected = [translator.features[i] for i in [1232]]
    # translator.features = selected

    # translator.log_feature_denotations(config.feature_denotation_filename, selected)

    translator.run()

    return dict(cnf_translator=translator)


def run_solver(config, data):
    # solution = solve(config.experiment_dir, config.cnf_filename, 'wpm3')
    # solution = solve(config.experiment_dir, config.cnf_filename, 'maxino')
    solution = solve(config.experiment_dir, config.cnf_filename, 'openwbo')
    if not solution.solved and solution.result == "UNSATISFIABLE":
        print_header("MAXSAT encoding is UNSATISFIABLE")
    else:
        print_header("MAXSAT solution with {} selected features found".format(solution.cost))
        data.cnf_translator.decode_solution(solution.assignment)

    return dict()


class Model(object):
    def __init__(self, cache):
        assert isinstance(cache, ExtensionCache)
        self.cache = cache

    def compute_feature_value(self, feature, state_id, substitution):
        try:
            return self.cache.feature_value(feature, state_id)
        except KeyError:
            value = feature.value(self.cache, state_id, substitution)
            self.cache.register_feature_value(feature, state_id, value)
            return value


class ModelTranslator(object):
    def __init__(self, features, states, goal_states, transitions, cache, cnf_filename):
        self.states = states
        self.goal_states = goal_states
        self.transitions = transitions
        self.state_ids = sorted(list(self.states.keys()))
        self.substitution = {}

        self.model, self.features = self.compute_feature_extensions(features, cache)
        self.writer = CNFWriter(cnf_filename)

        self.var_selected = None
        self.var_d1 = None
        self.var_d2 = None

        self.d1_distinguishing_features = dict()

        self.n_selected_clauses = 0
        self.n_d1_clauses = 0
        self.n_d2_clauses = 0
        self.n_bridge_clauses = 0
        self.n_goal_clauses = 0
        self.n_undistinguishable_state_pairs = 0

    def create_bridge_clauses(self, d1_literal, s, t):
        # If there are no transitions (t, t') in the sample set, then we do not post the bridge constraint.
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

    def compute_d1_distinguishing_features(self):
        # self.undistinguishable_state_pairs = []
        for s1, s2 in itertools.combinations(self.state_ids, 2):
            distinguishing = set()

            for f in self.features:
                x1 = self.model.compute_feature_value(f, s1, self.substitution)
                x2 = self.model.compute_feature_value(f, s2, self.substitution)
                if f.bool_value(x1) != f.bool_value(x2):  # f distinguishes s1 and s2
                    distinguishing.add(f)

            if not distinguishing:
                self.n_undistinguishable_state_pairs += 1
                # self.undistinguishable_state_pairs.append((s1, s2))

            self.d1_distinguishing_features[(s1, s2)] = distinguishing

    def run(self):
        print("Generating MAXSAT model from {} features".format(len(self.features)), flush=True)
        self.qchanges = qchanges = compute_qualitative_changes(self.transitions, self.features, self.model, self.substitution)

        self.compute_d1_distinguishing_features()

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
            for f in (x for x in d2_distinguishing_features if x not in self.d1_distinguishing_features[(s0, t0)]):
                forward_clause_literals.append(self.writer.literal(self.var_selected[f], True))
                self.writer.clause([d2_lit, self.writer.literal(self.var_selected[f], False)])
                self.n_d2_clauses += 1

            self.writer.clause(forward_clause_literals)
            self.n_d2_clauses += 1

        # Add the weighted clauses to minimize the number of selected features
        for feat_var in self.var_selected.values():
            self.writer.clause([self.writer.literal(feat_var, False)], weight=1)
            self.n_selected_clauses += 1

        self.report_stats()

        # self.debug_tests()

        self.writer.save()

        return self.writer.mapping

    def compute_feature_extensions(self, features, cache):
        """ Cache all feature denotations and prune those which have constant denotation at the same time """
        model = Model(cache)
        pruned = []
        for f in features:
            all_values = [model.compute_feature_value(f, s, self.substitution) for s in self.state_ids]
            if all_values.count(all_values[0]) != len(all_values):  # i.e. all elements are equal
                pruned.append(f)
            else:
                logging.debug("Feature \"{}\" has constant denotation ({}) over all states and will be ignored"
                              .format(f, all_values[0]))

        return model, pruned

    def decode_solution(self, assignment):
        varmapping = self.writer.mapping
        true_variables = set(varmapping[idx] for idx, value in assignment.items() if value is True)
        feature_mapping = {variable: feature for feature, variable in self.var_selected.items()}
        assert len(feature_mapping) == len(self.var_selected)
        selected_features = [feature_mapping[v] for v in true_variables if v in feature_mapping]
        print("Selected features: ")
        print('\n'.join(str(f) for f in selected_features))

    def report_stats(self):

        d1_distinguishing = self.d1_distinguishing_features.values()
        avg_num_d1_dist_features = sum(len(feats) for feats in d1_distinguishing) / len(d1_distinguishing)

        print_header("Max-sat encoding stats", 1)
        print_lines("Number of D1-undistinguishable state pairs: {}".format(self.n_undistinguishable_state_pairs), 1)
        print_lines("Avg. # of D1-distinguishing features: {:.2f}".format(avg_num_d1_dist_features), 1)
        print_lines("Clauses (possibly with repetitions):".format(), 1)
        print_lines("Selected: {}".format(self.n_selected_clauses), 2)
        print_lines("D1: {}".format(self.n_d1_clauses), 2)
        print_lines("D2: {}".format(self.n_d2_clauses), 2)
        print_lines("Bridge: {}".format(self.n_bridge_clauses), 2)
        print_lines("Goal: {}".format(self.n_goal_clauses), 2)
        print_lines("TOTAL: {}".format(self.n_selected_clauses + self.n_d1_clauses + self.n_d2_clauses +
                                       self.n_bridge_clauses + self.n_goal_clauses), 2)

    def log_features(self, feature_filename):
        logging.info("Pruned feature set saved in file {}".format(feature_filename))
        with open(feature_filename, 'w') as f:
            for i, feat in enumerate(self.features, 0):
                f.write("{}:\t{}\n".format(i, feat))
                # serialized = serialize_to_string(feat)
                # f.write("{}:\t{}\n".format(feat, serialized))
                # assert deserialize_from_string(serialized) == feat

    def log_feature_denotations(self, feature_denotation_filename, selected=None):
        selected = selected or self.features
        d = defaultdict(list)
        for (f, s), v in self.model.cache.feature_values.items():
            if f in selected:
                d[s].append((f, v))

        states = sorted(d.keys())
        with open(feature_denotation_filename, 'w') as file:
            for s in states:
                for (f, v) in sorted(d[s], key=lambda x: str(x[0])):
                    print("{} on state {}: {}".format(f, s, v), file=file)

    def debug_tests(self):
        undist = self.undistinguishable_state_pairs
        for s, t in undist:
            for s_prime in self.transitions[s]:
                some_d2_indistinguishable = False
                for t_prime in self.transitions[t]:
                    idx = compute_d2_index(s, s_prime, t, t_prime)
                    if idx not in self.d2_distinguished:
                        some_d2_indistinguishable = True

                if not some_d2_indistinguishable:
                    print("State pair {} violates bridge clauses".format((s, t)))
                    for t_prime in self.transitions[t]:
                        for f in self.features:
                            print("Delta[{}] for feature {}: {}".format((s, s_prime), f, self.qchanges[(s, s_prime, f)]))
                            print("Delta[{}] for feature {}: {}".format((t, t_prime), f, self.qchanges[(t, t_prime, f)]))
                    assert False


class Variable(object):
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name


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


class Clause(object):
    def __init__(self, literals, weight=math.inf):
        # assert all(isinstance(l, Literal) for l in literals)
        # self.literals = tuple(literals)
        self.literals = literals
        # assert len(set(literals)) == len(self.literals)  # Make sure all literals are unique
        self.weight = weight

    def __str__(self):
        return "{{{}}} [{}]".format(','.join(str(l) for l in self.literals), self.weight)

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
    w = "#TOP#" if weight is math.inf else weight
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

    def clause(self, literals, weight=math.inf):
        assert self.closed
        # self.clauses.add(Clause(literals=literals, weight=weight)) # Keeping the set in memory is expensive!
        self.num_clauses += 1
        self.accumulated_weight += weight if weight is not math.inf else 0

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
