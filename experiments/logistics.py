from sltp.util.misc import update_dict, extend_namer_to_all_features


def experiments():

    base = dict(
        domain_dir="logistics98",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
        feature_namer=feature_namer,
    )

    exps = dict()

    # Goal: arbitrary logistics goal
    exps["p1"] = update_dict(
        base,
        experiment_type='incremental',
        instances=['sample{}.pddl'.format(i) for i in [2]],
        test_instances=[
            "prob{:02d}.pddl".format(i) for i in range(2, 5)
        ],
        test_policy_instances=all_instances(),
        # num_states="until_first_goal",
        num_states="all",
        num_tested_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        max_concept_size=10,
        concept_generator=None,
        parameter_generator=None,
    )

    exps["p2"] = update_dict(exps["p1"], instances=['sample3.pddl'])

    exps["p2_max"] = update_dict(
        exps["p2"],

        goal_selector=goal_selector,
        create_goal_features_automatically=True,
    )


    exps["p1_p"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps


def feature_namer(feature):
    s = str(feature)
    base = {
        # "": "",
        "Exists(in,<universe>)": "num-loaded-packages",
        "And(Not(Equal(at_g,at)),obj)": "num-packages-not-in-destiny",
        "Exists(at_g,Exists(Inverse(at),airplane))": "num-packages-whose-destiny-has-an-airplane",
        "Exists(at_g,Exists(Inverse(at),<universe>))": "num-packages-with-destiny",
        "Exists(at_g,airport)": "num-packages-whose-destiny-has-an-airport",
        "Exists(in,airplane)": "num-packages-in-airplane",
        # "And(Equal(at_g,Inverse(at)),airport)": "",
        "And(Exists(at,airport),obj)": "num-packages-at-city-with-airport",
        # "Exists(at,Forall(Inverse(at),truck))": "",
    }

    return extend_namer_to_all_features(base).get(s, s)


def all_instances():
    return ["prob0{}.pddl".format(i) for i in range(1, 3)]


def goal_selector(lang):
    return "And(Equal(at_g,at),obj)"
