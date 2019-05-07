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
        test_policy_instances=all_test_instances(),
        num_states="until_first_goal",
        num_tested_states=20000,
        num_sampled_states=200,
        max_concept_size=8,
        distance_feature_max_complexity=8,
        cond_feature_max_complexity=8 + 2,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps


def all_test_instances():
    return ["problem{:02d}-full.pddl".format(i) for i in range(2, 12)]
