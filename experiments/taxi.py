#! /usr/bin/env python3
# -*- coding: utf-8 -*-
from sltp.util.misc import update_dict


def experiments():
    base = dict(
        domain_dir="taxi",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["simple"] = update_dict(
        base,
        instances=['instance_5.pddl', ],
        test_instances=["instance_7.pddl", ],
        num_states="until_first_goal",
        num_tested_states=50000,
        num_sampled_states=300,
        max_concept_size=6,
        distance_feature_max_complexity=6,
        feature_namer=feature_namer,
        feature_generator=expected_features,
        parameter_generator=None,
    )

    return exps


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


def expected_features(lang):
    return [
        "Bool[And(loct,locp)]",
        "Bool[And(locp,Nominal(inside_taxi))]",
        "Bool[And(locp,locp_g)]",
        "Dist[loct;adjacent;locp]",
        "Dist[loct;adjacent;locp_g]",
    ]
