#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from tarski.dl import PrimitiveConcept, ConceptCardinalityFeature, AndConcept, ExistsConcept, PrimitiveRole, StarRole, \
    InverseRole

from defaults import generate_experiment
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "gridworld"
    # domain_dir = "gripper-m"
    domain = "domain_strips.pddl"

    # This one overfits: in a 3x3 grid, with 2 booleans per dimension you can represent
    # the position
    sample_3 = dict(
        instances=["instance_strips_3.pddl"],
        complete_only_wrt_optimal=True,
        num_states=80, max_concept_size=6, max_concept_grammar_iterations=3,
        # distance_feature_max_complexity=6,
        concept_generator=None, parameter_generator=add_domain_parameters_strips)

    sample_5 = update_dict(sample_3, instances=["instance_strips_5.pddl"], num_states=100)
    sample_5_not_corner = update_dict(sample_5, instances=["instance_strips_5_not_corner.pddl"])
    # sample_5_not_corner = update_dict(sample_5, instances=["instance_strips_5_not_corner.pddl","instance_strips_5_not_corner_2.pddl"])

    s5_2inst = update_dict(sample_5, instances=["instance_strips_5.pddl", "instance_strips_5_not_corner.pddl"],
                           max_concept_size=10
                           # , feature_generator=generate_features_1
    )
    s5_with_feat = update_dict(sample_5, instances=["instance_strips_5.pddl"], feature_generator=generate_features_1)

    sample_5_rnd = update_dict(sample_5, num_states=1000, num_sampled_states=80, random_seed=12)
    sample_5_not_corner_rnd = update_dict(sample_5_not_corner, num_states=1000, num_sampled_states=80, random_seed=12)

    parameters = {
        "sample_3": sample_3,
        "sample_5": sample_5,
        "sample_5_not_corner": sample_5_not_corner,
        "sample_5_rnd": sample_5_rnd,
        "sample_5_not_corner_rnd": sample_5_not_corner_rnd,
        "s5_with_feat": s5_with_feat,
        "s5_2inst": s5_2inst,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def add_domain_parameters(language):
    # language.constant(2, "coordinate")  # x-goal coordinate
    # language.constant(3, "coordinate")  # x-goal coordinate
    # language.constant(10, "coordinate")  # grid limits!!
    # [language.constant(i, "coordinate") for i in range(1, 11)]
    return [language.constant(1, "coordinate")]  # grid limits!!


def add_domain_parameters_strips(language):
    return []


def generate_features_1(lang):
    obj_t = lang.Object
    maxpos = PrimitiveConcept(lang.get("maxpos"))
    xpos = PrimitiveConcept(lang.get("xpos"))
    xposg = PrimitiveConcept(lang.get("goal_xpos"))
    ypos = PrimitiveConcept(lang.get("ypos"))
    yposg = PrimitiveConcept(lang.get("goal_ypos"))
    succ = PrimitiveRole(lang.get("succ"))
    succ_i_s = StarRole(InverseRole(succ))
    succ_s = StarRole(succ)

    x_at_goal = AndConcept(xpos, xposg, "object")
    y_at_goal = AndConcept(ypos, yposg, "object")
    x_at_max = AndConcept(xpos, maxpos, "object")
    y_at_max = AndConcept(ypos, maxpos, "object")
    below_x = ExistsConcept(succ_s, xpos)
    below_xg = ExistsConcept(succ_s, xposg)
    below_y = ExistsConcept(succ_s, ypos)
    below_yg = ExistsConcept(succ_s, yposg)
    above_x = ExistsConcept(succ_i_s, xpos)
    above_y = ExistsConcept(succ_i_s, ypos)

    dist_x = AndConcept(below_x, below_xg, "object")
    dist_y = AndConcept(below_y, below_yg, "object")

    dist_x_true = AndConcept(above_x, below_xg, "object")
    dist_y_true = AndConcept(above_y, below_yg, "object")

    return [
        ConceptCardinalityFeature(x_at_goal),
        ConceptCardinalityFeature(y_at_goal),
        ConceptCardinalityFeature(x_at_max),
        ConceptCardinalityFeature(y_at_max),
        ConceptCardinalityFeature(below_x),
        ConceptCardinalityFeature(below_y),
        ConceptCardinalityFeature(above_x),
        ConceptCardinalityFeature(above_y),
        ConceptCardinalityFeature(dist_x),
        ConceptCardinalityFeature(dist_y),
        ConceptCardinalityFeature(dist_x_true),
        ConceptCardinalityFeature(dist_y_true),
    ]


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
