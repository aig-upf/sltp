from sltp.util.misc import update_dict
from sltp.util.names import blocksworld_names, blocksworld_parameters_for_clear, blocksworld_parameters_for_on


def experiments():
    base = dict(
        domain_dir="blocks",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=blocksworld_names,
    )

    exps = dict()

    exps["clear"] = update_dict(
        base,
        pipeline="transition_classifier",
        instances=["training_clear_5.pddl"],
        test_instances=[
            "instance_5_clear_x_1.pddl",
            "instance_5_clear_x_2.pddl",
        ],
        test_policy_instances=[
            "instance_5_clear_x_1.pddl",
            "instance_5_clear_x_2.pddl",
        ],
        num_states="all",
        max_concept_size=8,
        # concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        parameter_generator=blocksworld_parameters_for_clear,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete",
        use_equivalence_classes=True,
        use_feature_dominance=True,
    )

    exps["clear_fn"] = update_dict(
        exps["clear"],
        domain="domain_fstrips.pddl",
        test_domain="domain_fstrips.pddl",
        instances=["training_clear_5_fs.pddl"],
        test_instances=[],
        test_policy_instances=[],
    )

    exps["on"] = update_dict(
        base,
        pipeline="transition_classifier",
        # instances=["inst_on_x_y_16.pddl",
        #            "inst_on_x_y_14.pddl",
        #            "holding_a_b_unclear.pddl",
        #            ],
        instances=[
            "training_on_2.pddl",
            "training_on_3.pddl",
            "training_on_5.pddl",
        ],
        test_instances=["inst_on_x_y_7.pddl"],
        test_policy_instances=[
            "inst_on_x_y_7.pddl",
            "inst_on_x_y_14.pddl",
            "inst_on_x_y_16.pddl",
            "instance_9_on_x_y_1.pddl",
            "instance_3_on_x_y_2.pddl",
        ],
        # num_states=2000,
        # num_sampled_states=[50, 50, 1],
        # sampling="all",
        num_states="all",
        max_concept_size=8,
        # concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        # enforce_features=get_on_x_y_feature,
        # feature_generator=features_clear_x,
        parameter_generator=blocksworld_parameters_for_on,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete",
        use_equivalence_classes=True,
        use_feature_dominance=True,
    )

    exps["on_fn"] = update_dict(
        exps["on"],
        domain="domain_fstrips.pddl",
        test_domain="domain_fstrips.pddl",
        instances=["training_on_5_fs.pddl"],
        test_instances=[],
        test_policy_instances=[],
    )

    exps["all"] = update_dict(
        base,
        pipeline="transition_classifier",

        instances=[
            "probBLOCKS-5-0.pddl"
        ],
        test_instances=[

        ],
        test_policy_instances=[
            "probBLOCKS-6-0.pddl",
            "probBLOCKS-6-1.pddl",
            "probBLOCKS-6-2.pddl",
            "probBLOCKS-7-0.pddl",
            "probBLOCKS-7-1.pddl",
            "probBLOCKS-7-2.pddl",

        ],
        # num_states=2000,
        # num_sampled_states=[50, 50, 1],
        # sampling="all",
        num_states="all",
        max_concept_size=8,
        concept_generation_timeout=10,  # in seconds
        concept_generator=None,
        parameter_generator=None,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete",
        use_equivalence_classes=True,
    )

    return exps
