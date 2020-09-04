from sltp.util.misc import update_dict


def experiments():
    base = dict(
        domain_dir="visitall-opt11-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        # feature_namer=gripper_names,
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

    exps["small"] = update_dict(
        base,
        instances=[
            'problem03-full.pddl',
            # 'problem04-full.pddl',
            # 'problem05-full.pddl',
        ],
        test_instances=[],
        test_policy_instances=all_test_instances(),

        max_concept_size=8,
        distance_feature_max_complexity=8,
        cond_feature_max_complexity=8,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        # use_incremental_refinement=True,

        # feature_generator=debug_features,
    )

    return exps


def all_test_instances():
    return ["problem{:02d}-full.pddl".format(i) for i in range(2, 12)]


def debug_features(lang):
    nallbelow_wellplaced = "Num[Forall(Star(on),Equal(on_g,on))]"
    ontarget = "Num[Equal(on_g,on)]"
    nclear = f"Num[clear]"
    return [nallbelow_wellplaced, nclear, ontarget]
