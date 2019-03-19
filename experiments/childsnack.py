from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="childsnack-opt14-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=None,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            'sample01.pddl',
        ],
        test_instances=[
            # 'pfile01-001.pddl',
        ],
        num_states=200000,
        num_tested_states=20000,
        num_sampled_states=300,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps
