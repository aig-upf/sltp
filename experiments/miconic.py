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
        num_states=200000,
        num_tested_states=50000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        batch_refinement_size=5,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps
