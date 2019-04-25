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
        instances=[
            'sample01.pddl',
        ],
        test_instances=[
            'child-snack_pfile01.pddl',
        ],
        num_states=200000,
        num_tested_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=12, concept_bound_step=1,
        clean_workspace=False,
        concept_generator=None,
        parameter_generator=None,
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
        num_states=100,
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
        "served": "n-served",
        "notexist": "n-prepared-sndw",
        "no_gluten_sandwich": "n-no_gluten_sandwich",
        "at_kitchen_sandwich": "n-sdwch-at-kitchen",
        "Exists(ontray,<universe>)": "n-sdwch-ontray",
        "And(allergic_gluten,served)": "n-allergic-served",
        "Exists(at,Nominal(kitchen))": "n-trays-on-kitchen",
        "And(Not(served),not_allergic_gluten)": "n-normal-unserved",
    }

    return extend_namer_to_all_features(base).get(s, s)
