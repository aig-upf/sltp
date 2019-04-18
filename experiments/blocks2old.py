#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.incremental import IncrementalExperiment

from sltp.util.defaults import generate_experiment
from common import ijcai_paper_bw_feature_namer, \
    add_bw_domain_parameters_2
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "blocks"
    domain = "domain.pddl"

    exps = dict()

    exps["on_x_y"] = dict(
        instances=["inst_on_x_y_14.pddl"],
        num_states=40000, max_width=[-1],
        num_sampled_states=[500],
        # Note: Testing here is not simple, as we'd want to test only when X and Y are on different towers
        # test_domain=domain, test_instances=["inst_on_x_y_16.pddl"], num_tested_states=10000,
        complete_only_wrt_optimal=True,
        sampling="all",
        max_concept_size=8,
        concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,
    )

    exps["on_x_y_gc"] = update_dict(exps["on_x_y"], parameter_generator=None)

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
        instances=["probBLOCKS-8-0.pddl", ],
        # instances=["probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl",]
        test_instances=[
            "probBLOCKS-4-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=100000,
        num_sampled_states=None,  # To take as a basis all states
        initial_sample_size=100,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        batch_refinement_size=5,
        # quiet=True,
        clean_workspace=False,
    )

    exps["small_inc"] = update_dict(
        exps["arbitrary1"],
        experiment_class=IncrementalExperiment,
        instances=["probBLOCKS-4-0.pddl", ],
        test_instances=["probBLOCKS-3-0.pddl", ],
        num_states=10000,
        num_sampled_states=None,  # To take as a basis all states
        initial_sample_size=5,
        initial_concept_bound=5, max_concept_bound=12, concept_bound_step=1,
        batch_refinement_size=1,
        # quiet=True,
        clean_workspace=False,
    )

    exps["7blocks_inc"] = update_dict(
        exps["arbitrary1"],
        experiment_class=IncrementalExperiment,
        instances=["probBLOCKS-7-0.pddl", ],
        # instances=["probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl",]
        test_instances=[
            "probBLOCKS-5-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=50000,  # 50K is the minimum we'll need for probBLOCKS-7-0.pddl
        num_sampled_states=None,  # To take as a basis all states
        initial_sample_size=50,
        initial_concept_bound=5, max_concept_bound=14, concept_bound_step=1,
        batch_refinement_size=5,
        # quiet=True,
        clean_workspace=False,
    )

    exps["8blocks_inc"] = update_dict(
        exps["arbitrary1"],
        experiment_class=IncrementalExperiment,
        instances=["probBLOCKS-8-0.pddl", ],
        # instances=["probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl",]
        test_instances=[
            "probBLOCKS-4-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=600000,  # 600K is the minimum we'll need for probBLOCKS-8-0.pddl
        num_sampled_states=None,  # To take as a basis all states
        initial_sample_size=100,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        batch_refinement_size=5,
        # quiet=True,
        clean_workspace=False,
    )

    exps["7blocks_inc_k0_12"] = update_dict(
        exps["7blocks_inc"],
        initial_concept_bound=12, max_concept_bound=14, concept_bound_step=2,
    )

    exps["5blocks_inc_k0_12"] = update_dict(
        exps["7blocks_inc_k0_12"],
        initial_concept_bound=7, max_concept_bound=10, concept_bound_step=1,
        instances=["probBLOCKS-5-0.pddl", ],
        initial_sample_size=1000,
    )

    exps["5blocks_k7"] = dict(
        instances=["probBLOCKS-5-0.pddl"],
        test_domain=domain,
        test_instances=[
            "probBLOCKS-7-0.pddl",
        ],
        num_states=1000,
        num_tested_states=50000,
        num_sampled_states=None,
        complete_only_wrt_optimal=True,
        max_concept_size=7,
        feature_namer=ijcai_paper_bw_feature_namer,
    )

    exps["7blocks_manual"] = update_dict(
        exps["7blocks_inc"],
        feature_generator=genfeatures,
    )

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
    parameters["domain_dir"] = parameters.get("domain_dir", domain_dir)
    parameters["domain"] = parameters.get("domain", domain)
    return generate_experiment(**parameters)


def genfeatures():
    return ["Boolean[holding]",
            "Numerical[And(Forall(on_g,And(Equal(Star(on),Star(on_g)),clear)),clear)]",
            ]


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
