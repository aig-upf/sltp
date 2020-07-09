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
        # instances=["sample-2balls.pddl", "sample-small.pddl"],
        instances=["instance_5_clear_x_1.pddl"],
        test_instances=["instance_5_clear_x_2.pddl"],
        test_policy_instances=["instance_5_clear_x_2.pddl"],
        # num_states=2000,
        # num_sampled_states=300,
        num_states="all",
        max_concept_size=8,
        # concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        parameter_generator=blocksworld_parameters_for_clear,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete"
    )

    exps["on"] = update_dict(
        base,
        pipeline="transition_classifier",
        # instances=["inst_on_x_y_16.pddl",
        #            "inst_on_x_y_14.pddl",
        #            "holding_a_b_unclear.pddl",
        #            ],
        instances=["inst_on_x_y_5.pddl"],
        # We cannot test this, since works only for states where A and B are on diff towers
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
        optimal_selection_strategy="complete"
    )


    return exps
