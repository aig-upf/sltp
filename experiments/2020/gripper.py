from sltp.util.misc import update_dict
from sltp.util.names import gripper_names, gripper_parameters


def experiments():
    base = dict(
        domain_dir="gripper",
        domain="domain.pddl",
    )

    exps = dict()

    exps["small"] = update_dict(
        base,
        pipeline="transition_classifier",
        instances="sample-small.pddl",
        test_domain="domain.pddl",
        test_instances=["prob03.pddl", "prob04.pddl"],
        num_states=100,
        max_concept_size=10,
        max_concept_grammar_iterations=3,
        concept_generator=None,
        parameter_generator=gripper_parameters,
        feature_namer=gripper_names,
    )

    return exps
