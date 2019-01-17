#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "logistics98"
    domain = "domain.pddl"

    exps = dict()

    exps["sample1"] = dict(
        instances="sample1.pddl",
        num_states=100, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    exps["sample2"] = update_dict(
        exps["sample1"],
        instances="sample2.pddl",
        test_domain=domain, test_instances=["prob02.pddl"],
        pipeline="maxsat_poly",
    )
    exps["prob01"] = update_dict(exps["sample1"], instances="prob01.pddl")

    exps["prob01_blai"] = update_dict(
        exps["prob01"], pipeline="maxsat_poly",)

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
    return generate_experiment(domain_dir, domain, **parameters)


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
