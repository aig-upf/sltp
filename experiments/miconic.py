from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="miconic",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
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
        test_policy_instances=all_test_instances(),
        num_states="until_first_goal",
        num_tested_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        batch_refinement_size=5,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps


def all_test_instances():
    instances = []
    for i in range(1, 31, 3):  # jump 3-by-3 to have fewer instances
        for j in range(0, 5):  # Each x has 5 subproblems
            instances.append("s{}-{}.pddl".format(i, j))
    return instances
