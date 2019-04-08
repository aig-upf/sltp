#! /usr/bin/env python3
# -*- coding: utf-8 -*-
from common import ijcai_paper_bw_feature_namer, add_bw_domain_parameters_2
from sltp.util.misc import update_dict


def genfeatures():
    return [
        "Boolean[holding]",
        "Numerical[And(Forall(on_g,And(Equal(Star(on),Star(on_g)),clear)),clear)]",
    ]


def experiments():

    base = dict(
        domain_dir="blocks",
        domain="domain.pddl",
        # test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["on_x_y"] = update_dict(
        base,
        instances=["inst_on_x_y_14.pddl"],
        num_states=40000, max_width=[-1],
        num_sampled_states=[500],
        # Note: Testing here is not simple, as we'd want to test only when X and Y are on different towers
        # test_domain=domain, test_instances=["inst_on_x_y_16.pddl"], num_tested_states=10000,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,
    )

    exps["on_x_y_gc"] = update_dict(exps["on_x_y"], parameter_generator=None)

    # A small incremental experiment mostly for testing purposes
    exps["small_inc"] = update_dict(
        base,
        experiment_type='incremental',
        instances=["probBLOCKS-4-0.pddl", ],
        test_instances=["probBLOCKS-5-0.pddl", ],
        num_states=10000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=5,
        initial_concept_bound=5, max_concept_bound=12, concept_bound_step=1,
        batch_refinement_size=1,
        clean_workspace=False,
    )

    exps["arbitrary1"] = update_dict(
        base,
        instances=["probBLOCKS-4-0.pddl"],
        test_domain="domain.pddl",
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
        experiment_type='incremental',
        instances=["probBLOCKS-7-0.pddl", ],
        # instances=["probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl",]
        test_instances=[
            "probBLOCKS-4-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=50000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        batch_refinement_size=5,
        # quiet=True,
        clean_workspace=False,
    )

    return exps
