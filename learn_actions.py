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
import math
import sys
from signal import signal, SIGPIPE, SIG_DFL
from time import strftime, gmtime

from extensions import ExtensionCache
from utils import bootstrap
import features

signal(SIGPIPE, SIG_DFL)


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

    return qchanges


def create_variables(all_features, transitions, state_ids, writer):
    selected_variables = {feat: writer.variable("selected({})".format(feat)) for feat in all_features}

    d1_variables = {(s1, s2): writer.variable("D1[{}, {}]".format(s1, s2)) for s1, s2 in
                    itertools.combinations(state_ids, 2)}

    d2_variables = dict()

    for s, t in itertools.combinations(transitions, 2):
        for s_prime, t_prime in itertools.product(transitions[s], transitions[t]):
            idx = compute_d2_index(s, s_prime, t, t_prime)
            varname = "D2[({},{}), ({},{})]".format(*idx)
            d2_variables[idx] = writer.variable(varname)
    return selected_variables, d1_variables, d2_variables


def main(args):

    all_features, states, transitions, cache = features.main(args)

    # TODO REMOVE UNIVERSE AND EMPTY FEATURES

    model = Model(cache)
    substitution = {}
    qchanges = compute_qualitative_changes(transitions, all_features, model, substitution)

    state_ids = sorted(list(states.keys()))

    writer = CNFWriter()
    selected_variables, d1_variables, d2_variables = create_variables(all_features, transitions, state_ids, writer)

    # exp = sp.Cnf()

    for (s0, s1, t0, t1), d2_var in d2_variables.items():
        d2_distinguishing_features = []  # all features that d2-distinguish the current transition
        for f in all_features:
            if qchanges[(s0, s1, f)] != qchanges[(s0, s1, f)]:
                d2_distinguishing_features.append(f)

        # D2(s0,s1,t0,t2) iff OR_f selected(f), where f ranges over features that d2-distinguish the transition
        d2_lit = Literal(d2_var)
        forward_clause_literals = [-d2_lit]
        for f in d2_distinguishing_features:
            # big_or = big_or | selected_variables[f]
            forward_clause_literals.append(Literal(selected_variables[f]))
            writer.clause([d2_lit, -Literal(selected_variables[f])])

        writer.clause(forward_clause_literals)

    for s1, s2 in itertools.combinations(state_ids, 2):
        d1_variable = d1_variables[(s1, s2)]

        d1_distinguishing_features = []

        for f in all_features:
            x1 = model.compute_feature_value(f, s1, substitution)
            x2 = model.compute_feature_value(f, s2, substitution)
            if x1 != x2:  # f distinguishes s1 and s2
                d1_distinguishing_features.append(f)

        # Post the constraint: D1(si, sj) <=> OR active(f), where the OR ranges over all
        # those features f that tell apart si from sj
        d1_lit = Literal(d1_variable)
        forward_clause_literals = [-d1_lit]

        for f in d1_distinguishing_features:
            # big_or = big_or | selected_variables[f]
            forward_clause_literals.append(Literal(selected_variables[f]))
            writer.clause([d1_lit, -Literal(selected_variables[f])])

        writer.clause(forward_clause_literals)

        for s1_prime in transitions[s1]:  # will be empty set if not initialized, which is ok
            forward_clause_literals = [Literal(d1_variable)]

            for s2_prime in transitions[s2]:
                idx = compute_d2_index(s1, s1_prime, s2, s2_prime)
                forward_clause_literals.append(-Literal(d2_variables[idx]))

            writer.clause(forward_clause_literals)

        for s2_prime in transitions[s2]:  # will be empty set if not initialized, which is ok
            forward_clause_literals = [Literal(d1_variable)]
            for s1_prime in transitions[s1]:
                idx = compute_d2_index(s1, s1_prime, s2, s2_prime)
                forward_clause_literals.append(-Literal(d2_variables[idx]))
            writer.clause(forward_clause_literals)

        # Add the weighted clauses to minimize the number of selected features
        for feat_var in selected_variables.values():
            writer.clause([-Literal(feat_var)], weight=1)

    writer.save()


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
    def __init__(self):
        pass


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
        assert all(isinstance(l, Literal) for l in literals)
        self.literals = tuple(literals)
        assert len(set(literals)) == len(self.literals)  # Make sure all literals are unique
        self.weight = weight

    def __str__(self):
        return "{{{}}}".format(','.join(str(l) for l in self.literals))

    def cnf_line(self, top, variable_index):
        # w <literals> 0
        w = top if self.weight is math.inf else self.weight
        literals = " ".join(l.to_cnf(variable_index) for l in self.literals)
        return "{} {} 0".format(w, literals)

    def __eq__(self, other):
        return self.literals == other.literals

    def __hash__(self):
        return hash(self.literals)


class CNFWriter(object):
    def __init__(self):
        self.variables = dict()
        self.clauses = set()

    def variable(self, name):
        return self.variables.setdefault(name, Variable(name))

    def clause(self, literals, weight=math.inf):
        self.clauses.add(Clause(literals=literals, weight=weight))

    def save(self, filename="model.cnf"):
        numvars = len(self.variables)
        numclauses = len(self.clauses)
        top = sum(c.weight for c in self.clauses if c.weight is not math.inf) + 1
        print("Writing model to file \"{}\"".format(filename))
        print("Model has {} variables and {} clauses. Top weight is {}".format(numvars, numclauses, top))

        variable_index = {var: i for i, var in enumerate(self.variables.values(), start=1)}
        real_clause_printer = lambda c: c.cnf_line(top, variable_index)
        self._save(filename, numvars, numclauses, top, real_clause_printer)

        debug = True
        if debug:
            dfilename = "{}.txt".format(filename)
            debug_clause_printer = lambda c: str(c)
            self._save(dfilename, numvars, numclauses, top, debug_clause_printer)

    def _save(self, filename, numvars, numclauses, top, clause_printer):
        with open(filename, "w") as file:
            print("c WCNF model generated on {}".format(strftime("%Y%m%d %H:%M:%S", gmtime())), file=file)
            # p wcnf nbvar nbclauses top
            print("p wcnf {} {} {}".format(numvars, numclauses, top), file=file)
            for clause in self.clauses:
                print(clause_printer(clause), file=file)


if __name__ == "__main__":
    _args = bootstrap(sys.argv[1:])
    main(_args)


