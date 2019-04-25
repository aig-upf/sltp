from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="grid",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
        instances=[
            "prob01.pddl"
        ],
        test_instances=[
            'prob02.pddl',
        ],
        num_states=200000,
        num_tested_states=50000,
        num_sampled_states=None,  # Take all expanded states into account
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps
