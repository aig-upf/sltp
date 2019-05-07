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
        test_policy_instances=all_test_instances(),
        num_states="until_first_goal",
        num_tested_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=12, concept_bound_step=1,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps


def all_test_instances():
    return ["prob0{}.pddl".format(i) for i in range(1, 6)]
