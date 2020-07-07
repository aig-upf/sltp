from sltp.util.misc import update_dict
from sltp.util.names import gripper_names, gripper_parameters


def experiments():
    base = dict(
        # domain_dir="gripper-m",
        domain_dir="gripper",
        domain="domain.pddl",
    )

    exps = dict()

    exps["small"] = update_dict(
        base,
        pipeline="transition_classifier",
        # instances=["sample-2balls.pddl"],
        instances=["sample-small.pddl", "prob01.pddl"],
        test_domain="domain.pddl",
        test_instances=["prob03.pddl"],
        # test_instances=[],
        num_states="all",
        max_concept_size=10,
        concept_generator=None,
        parameter_generator=gripper_parameters,
        # parameter_generator=None
        feature_namer=gripper_names,
        maxsat_encoding="separation",
        complete_only_wrt_optimal=True,
        prune_redundant_states=False,
        optimal_selection_strategy="complete"
    )

    return exps
