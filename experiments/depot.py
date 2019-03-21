from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="depot",
        domain="domain.pddl",
        test_domain="domain.pddl",
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            'p01.pddl',
            # 'p02.pddl',
            # 'p03.pddl',
        ],
        test_instances=[
            # 'p03.pddl',
            # 'p04.pddl',
            # 'p05.pddl',
            # 'p06.pddl',
            # 'p07.pddl',
            'p08.pddl',
        ],
        num_states=2000,
        num_tested_states=20000,
        num_sampled_states=500,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p1_inc"] = update_dict(
        exps['p1'],
        instances=[
            'p01.pddl',
        ],
        experiment_type='incremental',
        num_states=10000,
        initial_sample_size=100,
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        batch_refinement_size=1,
    )

    exps["p1_inc_deb"] = update_dict(
        exps['p1_inc'],
        num_states=2000,
        num_tested_states=500,
        initial_sample_size=10,
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        batch_refinement_size=1,
    )

    return exps
