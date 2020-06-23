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
        instances=[
            "prob-6_4_10.pddl"
        ],
        test_instances=[
            # "prob-10-10-10-1540903568.pddl"
        ],
        num_states=20000,
        max_concept_size=12,
        concept_generation_timeout=60,  # in seconds
        concept_generator=None,
        parameter_generator=None,
        maxsat_encoding="separation",
    )

    return exps
