from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="depot",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
        instances=[
            'p01.pddl',
            'p02.pddl',
            # 'p03.pddl',
        ],
        test_instances=[
            # 'p03.pddl',
            # 'p04.pddl',
            # 'p05.pddl',
            # 'p06.pddl',
            'p07.pddl',
            'p08.pddl',
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
    return ["p{:02d}.pddl".format(i) for i in range(1, 23)]

