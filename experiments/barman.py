from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="barman-opt11-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=['sample01.pddl', ],
        test_instances=["pfile01-00{}.pddl".format(i) for i in (1, 2, 3)],
        num_states=200000,
        num_tested_states=50000,
        num_sampled_states=300,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_inc"] = update_dict(
        exps["p1"],
        experiment_type='incremental',
        test_instances=["pfile01-00{}.pddl".format(i) for i in (1, 2, 3)],
        num_states=200000,
        num_tested_states=50000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        clean_workspace=False,
    )

    exps["p1_p"] = update_dict(exps["p1"], pipeline="maxsat_poly")

    return exps
