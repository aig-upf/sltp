from sltp.util.misc import update_dict, extend_namer_to_all_features


def experiments():

    base = dict(
        domain_dir="childsnack-opt14-strips",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
        feature_namer=feature_namer,
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


def feature_namer(feature):
    s = str(feature)
    base = {
        # "": "",
        "served": "num-served-children",
        "And(Not(served),child)": "num-unserved-children",
        "notexist": "num-unprepared-sandwiches",
        "no_gluten_sandwich": "num-sandwiches-wo-gluten",
        "at_kitchen_bread": "num-breads-at-kitchen",
        "at_kitchen_content": "num-fillings-at-kitchen",
        "at_kitchen_sandwich": "num-sandwiches-at-kitchen",
        "Exists(ontray,<universe>)": "num-sandwiches-on-some-tray",
        "And(allergic_gluten,served)": "num-allergic-served",
        "Exists(at,Nominal(kitchen))": "num-trays-on-kitchen",
        "And(Not(served),not_allergic_gluten)": "num-unallergic-unserved",
        "And(Not(served),allergic_gluten)": "num-allergic-unserved",
        "And(Not(no_gluten_content),content-portion)": "num-gluten-free-fillings",
        "And(Not(no_gluten_bread),bread-portion)": "num-gluten-free-breads",
        "And(Not(no_gluten_sandwich),sandwich)": "num-sandwiches-with-gluten",
        "Exists(at,Exists(Inverse(waiting),<universe>))": "num-trays-on-place-with-some-child",
        "Exists(ontray,Exists(at,Nominal(kitchen)))": "num-sandwiches-on-some-tray-in-kitchen",
    }

    return extend_namer_to_all_features(base).get(s, s)


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
