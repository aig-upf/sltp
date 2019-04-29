#! /usr/bin/env python3
# -*- coding: utf-8 -*-
from sltp.util.misc import update_dict
from tarski.dl import PrimitiveConcept, GoalConcept, PrimitiveRole


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
        max_concept_size=8,
        distance_feature_max_complexity=8,
        feature_namer=feature_namer,
        concept_generator=build_expected_concepts,
        parameter_generator=None,
    )

    return exps


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


def build_expected_concepts(lang):
    """  """
    # obj_t = lang.Object
    loct, locp, adj = lang.get("loct", "locp", "adjacent")
    concepts = [PrimitiveConcept(loct), PrimitiveConcept(locp), GoalConcept(locp)]
    roles = [PrimitiveRole(adj)]
    return [], concepts, roles  # atoms, concepts, roles
