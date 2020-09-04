from sltp.util.misc import update_dict
from sltp.util.names import miconic_names


def experiments():
    base = dict(
        domain_dir="miconic",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=miconic_names,
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
            's2-0.pddl',
            # 's2-1.pddl',
            # 's2-2.pddl',
            # 's2-3.pddl',
            's3-0.pddl',
        ],
        test_instances=[
        ],
        test_policy_instances=all_test_instances(),

        max_concept_size=10,
        parameter_generator=None,
        use_equivalence_classes=True,
        # use_feature_dominance=True,
        use_incremental_refinement=True,
    )

    return exps


def all_test_instances():
    instances = []
    for i in range(1, 31, 3):  # jump 3-by-3 to have fewer instances
        for j in range(0, 5):  # Each x has 5 subproblems
            instances.append("s{}-{}.pddl".format(i, j))
    return instances
