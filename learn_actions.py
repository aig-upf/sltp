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
import sys
from signal import signal, SIGPIPE, SIG_DFL

import satispy as sp
from satispy.io import DimacsCnf

from extensions import ExtensionCache
from utils import bootstrap
import features

signal(SIGPIPE, SIG_DFL)


# class Variable(object):
#     def __init__(self, id_, name):
#         self.id = id_
#         self.name = name
#
#     def __str__(self):
#         return self.name


def cnf_tests():

    from satispy import Variable, Cnf
    from satispy.solver import Minisat

    v1 = Variable('v1')
    v2 = Variable('v2')
    v3 = Variable('v3')

    exp = v1 & v2 | v3

    io = DimacsCnf()
    res = io.tostring(exp)


    solver = Minisat()

    solution = solver.solve(exp)


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


def main(args):
    # cnf_tests()

    all_features, all_states, transitions, cache = features.main(args)
    active_variables = {feat: sp.Variable("act({})".format(feat)) for feat in all_features}

    state_ids = sorted(list(all_states.keys()))
    d1_variables = {(s1, s2): sp.Variable("D1({}, {})".format(s1, s2)) for s1, s2 in
                    itertools.combinations(state_ids, 2)}

    d2_variables = dict()

    for s, t in itertools.combinations(transitions, 2):
        for s_prime, t_prime in itertools.product(transitions[s], transitions[t]):
            idx = compute_d2_index(s, s_prime, t, t_prime)
            varname = "D2[({},{}), ({},{})]".format(*idx)
            d2_variables[idx] = sp.Variable(varname)

    # TODO REMOVE UNIVERSE AND EMPTY FEATURES


    model = Model(cache)
    substitution = {}

    qchanges = compute_qualitative_changes(transitions, all_features, model, substitution)

    exp = sp.Cnf()

    for (s0, s1, t0, t1), var in d2_variables.items():
        d2_distinguishing_features = []
        for f in all_features:
            if qchanges[(s0, s1, f)] != qchanges[(s0, s1, f)]:
                d2_distinguishing_features.append(f)

        big_or = sp.Cnf()
        for f in d2_distinguishing_features:
            big_or = big_or | active_variables[f]

        if len(d2_distinguishing_features) == 0:  # No distinguishing feature possible, hence D2 must be false
            exp = exp & -var
        else:
            exp = exp & (var >> big_or)
            exp = exp & (big_or >> var)

    for s1, s2 in itertools.combinations(state_ids, 2):
        d1_variable = d1_variables[(s1, s2)]

        d1_distinguishing_features = []

        for f in all_features:
            x1 = model.compute_feature_value(f, s1, substitution)
            x2 = model.compute_feature_value(f, s2, substitution)
            # if f distinguishes s1 and s2
            if x1 != x2:
                d1_distinguishing_features.append(f)

        # Post the constraint: D1(si, sj) <=> OR active(f), where the OR ranges over all
        # those features f that tell apart si from sj
        big_or = sp.Cnf()
        for f in d1_distinguishing_features:
            big_or = big_or | active_variables[f]
        if len(d1_distinguishing_features) == 0:  # No distinguishing feature possible, hence D1 must be false
            exp = exp & -d1_variable
        else:
            exp = exp & (d1_variable >> big_or)
            exp = exp & (big_or >> d1_variable)

        for s1_prime in transitions[s1]:  # will be empty set if not initialized, which is ok
            big_or = sp.Cnf()
            for s2_prime in transitions[s2]:
                idx = compute_d2_index(s1, s1_prime, s2, s2_prime)
                big_or = big_or | -d2_variables[idx]
            exp = exp & (-d1_variable >> big_or)

        for s2_prime in transitions[s2]:  # will be empty set if not initialized, which is ok
            big_or = sp.Cnf()
            for s1_prime in transitions[s1]:
                idx = compute_d2_index(s1, s1_prime, s2, s2_prime)
                big_or = big_or | -d2_variables[idx]
            exp = exp & (-d1_variable >> big_or)

        with open("model.cnf", "w") as f:
            io = DimacsCnf()
            f.write(io.tostring(exp))
            f.flush()


def store_model(clauses):
    with open("model.cnf", "w") as f:
        io = DimacsCnf()
        f.write(io.tostring(clauses))
        f.flush()

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


if __name__ == "__main__":
    _args = bootstrap(sys.argv[1:])
    main(_args)

