from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="satellite",
        domain="domain.pddl",
        test_domain="domain.pddl",
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            'p01-pfile1.pddl',
        ],
        test_instances=[
            'p05-pfile5.pddl',
        ],
        num_states=2000,
        num_tested_states=10000,
        num_sampled_states=200,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps
