from sltp.util.misc import update_dict
from sltp.util.names import logistics_names


def experiments():
    base = dict(
        # domain_dir="gripper-m",
        domain_dir="logistics98",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=logistics_names,
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

    # Goal: arbitrary logistics goal
    exps["small"] = update_dict(
        base,
        instances=['sample{}.pddl'.format(i) for i in [2]],
        # test_instances=["prob{:02d}.pddl".format(i) for i in range(2, 5)],
        test_instances=[],
        test_policy_instances=all_instances(),

        # num_states="until_first_goal",
        distance_feature_max_complexity=8,
        max_concept_size=8,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        use_incremental_refinement=True,
    )

    exps["p2_max"] = update_dict(
        exps["small"],

        goal_selector=goal_selector,
        create_goal_features_automatically=True,
    )

    return exps


def all_instances():
    return [f"prob0{i}.pddl" for i in range(1, 3)]


def goal_selector(lang):
    return "And(Equal(at_g,at),obj)"
