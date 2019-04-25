from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="logistics98",
        domain="domain.pddl",
        test_domain="domain.pddl",
    )

    exps = dict()

    # Goal: arbitrary logistics goal
    exps["p1"] = update_dict(
        base,
        instances=[
            'sample2.pddl',
        ],
        test_instances=["prob0{}.pddl".format(i) for i in (1, 2, 3, 4, 5)],
        num_states="until_first_goal",
        num_tested_states=50000,
        num_sampled_states=300,
        complete_only_wrt_optimal=True,
        max_concept_size=10,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    exps["p2"] = update_dict(
        exps["p1"],
        instances="sample2.pddl",
        test_instances=["prob02.pddl"],)

    exps["p2_p"] = update_dict(
        exps["p2"], pipeline="maxsat_poly")

    return exps
