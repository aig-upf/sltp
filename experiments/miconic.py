from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="miconic",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=None,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            's2-0.pddl',
            # 's2-1.pddl',
            # 's2-2.pddl',
            # 's2-3.pddl',
            's3-0.pddl',
        ],
        test_instances=[
            's2-4.pddl',
            's3-1.pddl',
        ],
        num_states=2000,
        num_tested_states=20000,
        num_sampled_states=200,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps
