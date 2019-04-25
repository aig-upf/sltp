from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="visitall-opt11-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            'problem03-full.pddl',
        ],
        test_instances=[
            'problem04-full.pddl',
        ],
        num_states=2000,
        num_tested_states=20000,
        num_sampled_states=200,
        max_concept_size=8,
        distance_feature_max_complexity=8,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps
