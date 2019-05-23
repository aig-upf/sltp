from sltp.util.misc import update_dict, extend_namer_to_all_features


def experiments():

    base = dict(
        domain_dir="spanner-small",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
        feature_namer=feature_namer,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        # experiment_type='incremental',
        instances=[
            # "prob-4-4-3-1540907456.pddl",
            # "prob-6_4_2.pddl",
            "prob-6_4_10.pddl"
        ],
        test_instances=[
            # "prob-10-10-10-1540903568.pddl"
        ],
        test_policy_instances=all_test_instances(),
        # num_states="all",
        num_states=20000,
        num_sampled_states=None,  # Take all expanded states into account
        num_tested_states=20000,
        initial_sample_size=100, batch_refinement_size=5,
        initial_concept_bound=8, max_concept_bound=16, concept_bound_step=1,
        distance_feature_max_complexity=8,
        cond_feature_max_complexity=8 + 2,
        concept_generator=None,
        # goal_selector=goal_selector,
        # create_goal_features_automatically=True,

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


def feature_namer(feature):
    s = str(feature)
    base = {
        # "": "",
        "And(tightened_g,Not(tightened))": "n-untightened-nuts",
        "Exists(carrying,<universe>)": "n-carried-spanners",
        "Forall(Inverse(link),<empty>)": "first-cell",  # Neat!
        "Exists(at,Forall(Inverse(link),<empty>))": "n-things-on-first-cell",
        "And(Exists(at,Exists(Inverse(at),man)),Not(man))": "n-spanners-in-same-cell-as-man",
        "And(Exists(at,Exists(Inverse(at),man)),spanner)":  "n-spanners-in-same-cell-as-man",
        "Exists(at,Exists(link,Exists(Inverse(at),<universe>)))": "",
        "loose": "n-untightened-nuts",
        "Exists(at,Exists(link,Exists(Inverse(at),man)))": "n-spanners-on-cell-left-to-man",
        "": "",
        "": "",
    }
    return extend_namer_to_all_features(base).get(s, s)

# 	3. Num[Exists(at,Exists(link,Exists(Inverse(at),<universe>)))] [k=7, id=27] ???

# 	4. Bool[Exists(at,Exists(link,Forall(link,Forall(link,<empty>))))] [k=8, id=55]
# 	5. Num[Exists(at,Exists(link,Exists(Inverse(at),man)))] [k=8, id=63]
# 	6. Bool[And(Exists(at,Exists(link,Forall(Inverse(at),<empty>))),man)] [k=9, id=108]

# 	8. Bool[And(Exists(link,Forall(Inverse(at),man)),Forall(Inverse(link),<empty>))] [k=10, id=264]

# 	10. Num[Exists(at,Exists(link,Exists(link,Exists(Inverse(at),man))))] [k=10, id=280]


# 	3. Num[Not(Forall(at,Forall(link,Forall(Inverse(at),spanner))))] [k=9, id=115]
# 	4. Bool[And(Exists(at,Forall(link,Exists(Inverse(at),<universe>))),man)] [k=9, id=139]
# 	5. Bool[And(Exists(at,Exists(link,Forall(link,Exists(link,<universe>)))),man)] [k=10, id=297]

# 	6. Bool[And(Exists(at,Exists(link,Exists(Inverse(at),spanner))),spanner)] [k=10, id=308]
# 	9. Num[    Exists(at,Exists(link,Exists(link,Exists(Inverse(at),spanner))))] [k=10, id=274]

# 	7. Num[Exists(at,Exists(link,Exists(link,Exists(Inverse(at),spanner))))] [k=10, id=427]


def goal_selector(lang):
    # return "And(nut,tightened)"
    return "Not(loose)"


def all_test_instances():
    instances = []
    total = 1
    for i in range(1, 7):
        for _ in range(5):  # Each x has 5 subproblems
            instances.append("pfile0{}-{:03d}.pddl".format(i, total))
            total += 1
    assert total == 31
    return instances
