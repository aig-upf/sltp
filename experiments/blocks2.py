#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.incremental import IncrementalExperiment

from defaults import generate_experiment
from common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, ijcai_paper_bw_feature_namer, \
    add_bw_domain_parameters_2, build_on_x_y_feature_set, generate_features_n_ab, get_on_x_y_feature, \
    features_clear_x
from sltp.util.misc import update_dict
from sltp.util.serialization import deserialize_from_string


def experiment(experiment_name=None):
    domain_dir = "blocks"
    domain = "domain.pddl"

    exps = dict()

    exps["arbitrary1"] = dict(
        instances=["probBLOCKS-4-0.pddl"],
        test_domain=domain,
        test_instances=[
            "probBLOCKS-5-0.pddl",
            "probBLOCKS-5-1.pddl",
            "probBLOCKS-6-0.pddl",
            "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-0.pddl",
            "probBLOCKS-10-1.pddl",
            # "probBLOCKS-10-2.pddl",
            "probBLOCKS-14-0.pddl",
            # "probBLOCKS-14-1.pddl",
        ],
        num_states=2000,
        num_tested_states=50000,
        num_sampled_states=300,
        complete_only_wrt_optimal=True,
        max_concept_size=6,
        concept_generator=None,
        parameter_generator=None,
        feature_namer=ijcai_paper_bw_feature_namer,
    )

    exps["arbitrary1_inc"] = update_dict(
        exps["arbitrary1"],
        experiment_class=IncrementalExperiment,
        instances=["probBLOCKS-6-0.pddl", ],
        # instances=["probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl",]
        test_instances=[
            "probBLOCKS-4-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=10000,
        initial_sample_size=100,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        batch_refinement_size=5,
        # quiet=True,
        clean_workspace=False,
    )

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
    parameters["domain_dir"] = parameters.get("domain_dir", domain_dir)
    parameters["domain"] = parameters.get("domain", domain)
    return generate_experiment(**parameters)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
