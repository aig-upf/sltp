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
        instances=["prob-3-2-2.pddl"],
        test_instances=[
            "prob-4-3-3-1540907466.pddl"
        ],
        test_policy_instances=[
            "prob-10-10-10-1540903568.pddl",
            "prob-15-10-8-1540913795.pddl"
        ],
        num_states="all",
        max_concept_size=8,
        # concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        parameter_generator=None,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete",
        transition_classification_policy=debug_policy
    )

    return exps


def debug_policy():
    n_carried = "Num[Exists(Inverse(carrying),<universe>)]"  # K=3
    n_unreachable_locs = "Num[Exists(Star(link),Exists(Inverse(at),man))]"  # K=7
    n_tightened = "Num[tightened]"  # K=1
    n_spanners_same_loc_as_bob = "Num[And(Exists(at,Exists(Inverse(at),man)),spanner)]"  # K=8

    return [
        # picking a spanner when possible is always good:
        [(n_carried, 'INC')],

        # Tightening a nut when possible is always good:
        [(n_tightened, 'INC')],

        # Moving to the right when there's no more spanners in same location as bob is good:
        [(n_spanners_same_loc_as_bob, "=0"), (n_unreachable_locs, 'INC')],
    ]
