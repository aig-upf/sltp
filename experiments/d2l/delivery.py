
from sltp.util.misc import update_dict
from sltp.util.names import taxi_names


def experiments():
    base = dict(
        domain_dir="delivery",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=taxi_names,
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

        distinguish_goals=True,
    )

    exps = dict()

    exps["small"] = update_dict(
        base,
        instances=['instance_4_0.pddl'],
        test_instances=[],

        # Cannot verify, not in STRIPS
        # test_policy_instances=[f"instance_7_{i}.pddl" for i in range(0, 3)],

        max_concept_size=8,
        distance_feature_max_complexity=8,
        # feature_generator=expected_features,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        use_incremental_refinement=True,

        print_denotations=True,
    )

    return exps


def expected_features_wo_conditionals(lang):
    return [
        "Bool[And(loct,locp)]",
        "Bool[And(locp,Nominal(inside_taxi))]",
        "Bool[And(locp,locp_g)]",
        "Dist[loct;adjacent;locp]",
        "Dist[loct;adjacent;locp_g]",
    ]


def expected_features(lang):
    return [
        "Bool[And(locp,locp_g)]",  # Goal-distinguishing
        "Dist[loct;adjacent;locp]",  # Distance between taxi and passenger
        "Bool[And(loct,locp)]",
        "Bool[And(locp,Nominal(inside_taxi))]",
        "If{Bool[And(locp,Nominal(inside_taxi))]}{Dist[locp_g;adjacent;loct]}{Infty}",
    ]
