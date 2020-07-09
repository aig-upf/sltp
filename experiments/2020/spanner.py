from sltp.util.misc import update_dict
from sltp.util.names import spanner_names


def experiments():
    base = dict(
        domain_dir="spanner-small",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
        feature_namer=spanner_names,
    )

    exps = dict()

    exps["small"] = update_dict(
        base,
        pipeline="transition_classifier",
        instances=["prob-3-2-2-.pddl"],
        test_instances=[
            # "prob-10-10-10-1540903568.pddl"
        ],
        test_policy_instances=[
            # "prob-10-10-10-1540903568.pddl"
        ],
        num_states="all",
        max_concept_size=10,
        distance_feature_max_complexity=8,
        concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        parameter_generator=None,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete"
    )

    return exps
