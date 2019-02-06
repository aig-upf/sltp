#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "childsnack-opt14-strips"
    domain = "domain.pddl"

    exps = dict()

    exps["prob01"] = dict(
        instances="sample01.pddl",
        # test_domain=domain, test_instances=["pfile01-001.pddl"],
        num_states=200000, max_concept_size=8, max_concept_grammar_iterations=3,
        num_sampled_states=300,
        # distance_feature_max_complexity=8,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    exps["prob01_tg"] = update_dict(
        exps["prob01"], complete_only_wrt_optimal=True)

    exps["prob01_p"] = update_dict(
        exps["prob01"], pipeline="maxsat_poly")

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
