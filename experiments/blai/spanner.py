#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.util.defaults import generate_experiment
from sltp.util.misc import update_dict

#
# good abstraction:
#   1. #loose nuts: loose(?n)
#   2. #usable-spanners-being-carried: usable(?s) AND carry(m,?s)
#   3. #usable-spanners-not-being-carried: usable(?s) AND -carry(m,?s)
#   4. distance-to-usable-spanner-not-being-carried: dist(C1, R:C, C2), 
#   5. distance-to-gate: ??
#
# move-towards-usable-spanner: f3 > 0, f4 > 0; dec(f4)
# pick-spanner: f3 > 0, f4 = 0; dec(f3), inc(f2), inc*(f4)
# move-towards-gate: f5 > 0; dec(f5)
# tighten-nut: f1 > 0, f2 > 0; dec(f1), dec(f2)
#
# problem: abstraction may not be solvable by non-det
#          effects on f1 and f2 by tighted-nut
#

def experiment(experiment_name=None):
    domain_dir = "spanner-small"
    domain = "domain.pddl"

    exps = dict()

    exps["exp1"] = dict(
        instances="prob-4-4-3-1540907456.pddl",
        #instances="prob-3-3-3-1540903410.pddl",
        test_domain=domain,
        #test_instances=["prob-10-10-10-1540903568.pddl"],
        test_instances=[],
        num_states=523,
        max_concept_size=8,
        max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        concept_generator=None,
        parameter_generator=add_domain_parameters,  # No goal concepts!
        feature_namer=feature_namer,)

    exps["exp1_tg"] = update_dict(
        exps["exp1"],
        complete_only_wrt_optimal=True)

    exps["exp1_p"] = update_dict(
        exps["exp1"],
        pipeline="maxsat_poly")

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    return generate_experiment(domain_dir, domain, **exps[experiment_name])


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


def add_domain_parameters(language):
    return []


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
