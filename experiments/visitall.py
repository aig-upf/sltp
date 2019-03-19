#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.incremental import IncrementalExperiment

from defaults import generate_experiment
from common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, ijcai_paper_bw_feature_namer, \
    add_bw_domain_parameters_2, build_on_x_y_feature_set, generate_features_n_ab, get_on_x_y_feature, \
    features_clear_x
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "visitall-opt11-strips"
    domain = "domain.pddl"

    exps = dict()

    exps["p1"] = dict(
        instances=[
            'problem03-full.pddl',
        ],
        test_domain=domain,
        test_instances=[
            'problem04-full.pddl',
        ],
        num_states=2000,
        num_tested_states=20000,
        num_sampled_states=200,
        max_concept_size=8,
        distance_feature_max_complexity=8,
        complete_only_wrt_optimal=True,
        parameter_generator=None,
    )

    exps["p1_poly"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
    parameters["domain_dir"] = parameters.get("domain_dir", domain_dir)
    parameters["domain"] = parameters.get("domain", domain)
    return generate_experiment(**parameters)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
