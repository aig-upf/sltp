from sltp.util.misc import update_dict, extend_namer_to_all_features
from sltp.util.names import childsnack_names


def experiments():

    base = dict(
        domain_dir="childsnack-opt14-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
        feature_namer=childsnack_names,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
        # instances=['sample{:02d}.pddl'.format(i) for i in range(1, 5)],
        instances=['child-snack_pfile01-2.pddl'],
        test_instances=[
            # 'child-snack_pfile01-2.pddl',
        ],
        test_policy_instances=all_test_instances(),
        # num_states="until_first_goal",
        # num_states="all",
        num_states=20000,
        num_tested_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=12, concept_bound_step=1,
        clean_workspace=False,
        concept_generator=None,
        parameter_generator=None,

        goal_selector=goal_selector,
        create_goal_features_automatically=True,
    )

    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    exps["mini"] = update_dict(
        base,
        experiment_type='incremental',
        instances=[
            'sample_mini.pddl',
        ],
        test_instances=[
            'child-snack_pfile01.pddl',
        ],
        num_states="until_first_goal",
        num_tested_states=500,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=20, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        clean_workspace=False,
        concept_generator=None,
        parameter_generator=None,
    )

    return exps


def all_test_instances():
    return ['child-snack_pfile01.pddl', 'child-snack_pfile01-2.pddl', 'child-snack_pfile02.pddl',
            'child-snack_pfile04-2.pddl', 'child-snack_pfile05.pddl', 'child-snack_pfile07-2.pddl',
            'child-snack_pfile08.pddl', 'child-snack_pfile10-2.pddl', 'child-snack_pfile01.pddl',
            'child-snack_pfile03-2.pddl', 'child-snack_pfile04.pddl', 'child-snack_pfile06-2.pddl',
            'child-snack_pfile07.pddl', 'child-snack_pfile09-2.pddl', 'child-snack_pfile10.pddl',
            'child-snack_pfile02-2.pddl', 'child-snack_pfile03.pddl', 'child-snack_pfile05-2.pddl',
            'child-snack_pfile06.pddl', 'child-snack_pfile08-2.pddl', 'child-snack_pfile09.pddl']


def goal_selector(lang):
    return "served"
