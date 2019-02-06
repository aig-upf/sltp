#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "visitall-opt11"
    domain = "domain.pddl"

    exps = dict()

    exps["problem03full"] = dict(
        instances="problem03-full.pddl",
        test_domain=domain, test_instances=["problem04-full.pddl"],
        num_states=500, max_concept_size=8, max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    exps["problem03full_tg"] = update_dict(
        exps["problem03full"], complete_only_wrt_optimal=True)

    exps["problem03full_p"] = update_dict(
        exps["problem03full"], pipeline="maxsat_poly")

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    return generate_experiment(domain_dir, domain, **exps[experiment_name])


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
