from sltp.util.misc import update_dict
from sltp.util.names import blocksworld_names, blocksworld_parameters_for_clear, blocksworld_parameters_for_on


def experiments():
    base = dict(
        domain_dir="blocks",
        feature_namer=blocksworld_names,
        pipeline="transition_classifier",
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete",
        num_states="all",
        concept_generator=None,
        parameter_generator=None,
        v_slack=2,

        # concept_generation_timeout=120,  # in seconds
        maxsat_timeout=None,
    )

    exps = dict()

    strips_base = update_dict(
        base,
        domain="domain.pddl",
        test_domain="domain.pddl",
    )

    fn_base = update_dict(
        base,
        domain="domain_fstrips.pddl",
        test_domain="domain_fstrips.pddl",
    )

    strips_atomic_base = update_dict(
        base,
        domain="domain_atomic_move.pddl",
        test_domain="domain_atomic_move.pddl",
    )

    exps["clear"] = update_dict(
        strips_base,
        instances=["training_clear_5.pddl"],
        test_instances=[
            "instance_5_clear_x_1.pddl",
            # "instance_5_clear_x_2.pddl",
        ],
        test_policy_instances=[
            "testing_clear_10_0.pddl",
            "testing_clear_10_1.pddl",
        ],

        max_concept_size=8,
        parameter_generator=blocksworld_parameters_for_clear,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        use_incremental_refinement=True,
    )

    exps["clear_fn"] = update_dict(
        fn_base,
        instances=["training_clear_5_fs.pddl"],
        test_instances=[],
        test_policy_instances=[],
        max_concept_size=8,
        parameter_generator=blocksworld_parameters_for_clear,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        use_incremental_refinement=True,
    )

    exps["on"] = update_dict(
        strips_base,
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

        max_concept_size=8,
        # enforce_features=get_on_x_y_feature,
        parameter_generator=blocksworld_parameters_for_on,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        # use_incremental_refinement=True,
    )

    exps["on_fn"] = update_dict(
        fn_base,
        instances=["training_on_5_fs.pddl"],
        test_instances=[],
        test_policy_instances=[],

        max_concept_size=8,
        parameter_generator=blocksworld_parameters_for_on,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        # use_incremental_refinement=True,
    )

    exps["all_5"] = update_dict(
        strips_base,

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

        max_concept_size=10,
        use_incremental_refinement=True,
        use_equivalence_classes=True,
        use_feature_dominance=False,
    )

    exps["all_fn"] = update_dict(
        fn_base,
        instances=[
            # "training_singletower_5_fs.pddl",
            # "training_singletower_6_fs.pddl",
            # "training_arbitrary_5_fs.pddl",
            "training_arbitrary_6_fs.pddl",
            # "training_singletower_8_fs.pddl",
            # "training_singletower_10_fs.pddl",
        ],
        test_instances=[],
        test_policy_instances=[],

        # feature_generator=fs_debug_features,
        # transition_classification_policy=fs_debug_policy
        max_concept_size=8,
        use_equivalence_classes=True,
        use_incremental_refinement=True,
    )

    exps["all_fn_5"] = update_dict(
        exps["all_fn"],
        instances=[
            "training_arbitrary_5_fs.pddl",
            # "training_singletower_5_fs.pddl",
        ],
        max_concept_size=10,
        # feature_generator=fs_debug_features,
        use_incremental_refinement=False,
        use_equivalence_classes=True,
        use_feature_dominance=False,
    )

    exps["all_at_5"] = update_dict(
        strips_atomic_base,
        instances=[
            "training_arbitrary_5_atomic.pddl",
        ],
        test_instances=[],
        test_policy_instances=[
            "training_arbitrary_5_atomic.pddl",
            "testing_arbitrary_10_atomic.pddl",
            "testing_arbitrary_10_1_atomic.pddl",
            "testing_arbitrary_17-0_atomic.pddl",
            "testing_arbitrary_17-1_atomic.pddl",
        ],

        max_concept_size=10,
        # feature_generator=debug_features_at,
        use_incremental_refinement=True,
        use_equivalence_classes=True,
        use_feature_dominance=False,

        # force_zeros=True,
    )

    exps["all_at_testing"] = update_dict(
        strips_atomic_base,
        instances=[
            "training_arbitrary_5_atomic.pddl",
        ],
        test_instances=[],
        test_policy_instances=[
            "training_arbitrary_5_atomic.pddl",
            "testing_arbitrary_10_atomic.pddl",
            "testing_arbitrary_10_1_atomic.pddl",
            "testing_arbitrary_17-0_atomic.pddl",
            "testing_arbitrary_17-1_atomic.pddl",
        ],
        feature_generator=debug_features_at2,
        use_incremental_refinement=False,
        use_equivalence_classes=True,
        use_feature_dominance=False,
    )

    return exps


def fs_debug_policy():
    on = "loc"
    eqons = f'Equal({on}_g,{on})'
    wellplaced = f"Num[And({eqons},Forall(Star({on}),{eqons}))]"
    nclear = f"Num[clear]"

    return [
        # Increasing the number of well-placed blocks is always good
        [(wellplaced, 'INC')],

        # Increasing the # of clear blocks (i.e. by putting some block on the table) is always good as long
        # as that doesn't affect the number of well placed blocks
        [(wellplaced, 'NIL'), (nclear, 'INC')],
    ]


def fs_debug_features(lang):
    on = "loc"
    eqons = f'Equal({on}_g,{on})'
    wellplaced = f"Num[And({eqons},Forall(Star({on}),{eqons}))]"
    nclear = f"Num[clear]"
    return [wellplaced, nclear]


def fs_debug_features2(lang):
    nclear = "Num[clear]"
    sup_wp = "Num[Equal(Star(loc_g),Star(loc))]"
    ontarget = "Num[Equal(loc_g,loc)]"
    return [sup_wp, nclear, ontarget]


def debug_features_at(lang):
    # This alone is UNSAT
    on = "on"
    eqons = f'Equal({on}_g,{on})'
    wellplaced = f"Num[And({eqons},Forall(Star({on}),{eqons}))]"
    nclear = f"Num[clear]"
    return [wellplaced, nclear]


def debug_features_at2(lang):
    nallbelow_wellplaced = "Num[Forall(Star(on),Equal(on_g,on))]"
    ontarget = "Num[Equal(on_g,on)]"
    nclear = f"Num[clear]"
    return [nallbelow_wellplaced, nclear, ontarget]
