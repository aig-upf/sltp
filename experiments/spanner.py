from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="spanner-small",
        domain="domain.pddl",
        test_domain="domain.pddl",
        feature_namer=None,
    )

    exps = dict()

    exps["p1"] = update_dict(
        base,
        instances=[
            "prob-4-4-3-1540907456.pddl"
        ],
        test_instances=[
            "prob-10-10-10-1540903568.pddl"
        ],
        num_states=2000,
        num_tested_states=50000,
        num_sampled_states=500,
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        distance_feature_max_complexity=8,
        concept_generator=None,
        parameter_generator=add_domain_parameters,  # No goal concepts!
    )

    exps["p1"] = update_dict(
        exps["p1"], pipeline="maxsat_poly")

    return exps


def add_domain_parameters(language):
    return []
