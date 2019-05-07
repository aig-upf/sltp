from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="spanner-small",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
        instances=[
            "prob-4-4-3-1540907456.pddl",
            "prob-6_4_2.pddl",
            "prob-6_4_10.pddl"
        ],
        test_instances=[
            "prob-10-10-10-1540903568.pddl"
        ],
        test_policy_instances=all_test_instances(),
        num_states="until_first_goal",
        num_sampled_states=None,  # Take all expanded states into account
        num_tested_states=20000,
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        batch_refinement_size=5,
        distance_feature_max_complexity=8,
        concept_generator=None,

        # No goal concepts would be necessary here, but we let the spanner experiments run with them so that they are
        # homogeneous with the rest of experiments
        # parameter_generator=add_domain_parameters,  # This would prevent goal concepts from being generated
        parameter_generator=None,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps


def add_domain_parameters(language):
    return []


def all_test_instances():
    instances = []
    total = 1
    for i in range(1, 7):
        for _ in range(5):  # Each x has 5 subproblems
            instances.append("pfile0{}-{:03d}.pddl".format(i, total))
            total += 1
    assert total == 31
    return instances
