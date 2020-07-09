from sltp.util.misc import update_dict
from sltp.util.names import gripper_names, gripper_parameters


def experiments():
    base = dict(
        # domain_dir="gripper-m",
        domain_dir="gripper",
        domain="domain.pddl",
        feature_namer=gripper_names,
    )

    exps = dict()

    exps["small"] = update_dict(
        base,
        pipeline="transition_classifier",
        # instances=["sample-2balls.pddl", "sample-small.pddl"],
        instances=["prob01.pddl"],
        test_domain="domain.pddl",
        # test_instances=[],
        test_instances=[f"prob{i:02d}.pddl" for i in range(3, 11)],
        test_policy_instances=[f"prob{i:02d}.pddl" for i in range(3, 21)],
        num_states="all",
        max_concept_size=10,
        # concept_generation_timeout=120,  # in seconds
        concept_generator=None,
        parameter_generator=gripper_parameters,
        # parameter_generator=None
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete"
    )

    exps["quick"] = update_dict(
        exps["small"],
        test_instances=[],
        test_policy_instances=["prob04.pddl"],
    )

    return exps
