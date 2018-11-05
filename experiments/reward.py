#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from abstractions_defaults import generate_experiment
from common import update_dict
from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole, UniversalConcept, GoalConcept, MinDistanceFeature, ConceptCardinalityFeature, RestrictRole


def add_domain_parameters(language):
    return []


def build_expected_concepts(lang):
    obj_t = lang.Object

    top = UniversalConcept("object")

    current_cell = PrimitiveConcept(lang.get("at"))
    reward_cells = PrimitiveConcept(lang.get("reward"))
    blocked_cells = PrimitiveConcept(lang.get("blocked"))
    unblocked_cells = NotConcept(blocked_cells, obj_t)
    # visited_cells = PrimitiveConcept(lang.get("visited"))
    # unvisited_cells = NotConcept(visited_cells, obj_t)
    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    # unvisited_cells_with_reward = AndConcept(unvisited_cells, reward_cells, "cell")

    at_cell_with_reward = AndConcept(current_cell, reward_cells, "cell")  # X

    # concepts = [current_cell, unvisited_cells, unvisited_cells_with_reward, unblocked_cells]
    # concepts = [current_cell, unblocked_cells, reward_cells]
    # concepts = [current_cell, unblocked_cells, reward_cells, at_cell_with_reward]
    concepts = [top, current_cell, reward_cells, at_cell_with_reward]

    roles = [adjacent_role]

    return [], concepts, roles  # atoms, concepts, roles


def experiment(experiment_name=None):
    domain_dir = "grid-circles"
    domain = "domain.pddl"

    # One reason for overfitting: in a 3x3 grid, with 2 booleans per dimension you can perfectly represent any position
    sample_1x3 = dict(
        instances=["sample_1x3.pddl"],
        complete_only_wrt_optimal=True,
        num_states=1500, max_concept_size=8, max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        feature_namer=feature_namer,
        # relax_numeric_increase=True,
        # feature_generator=build_expected_features,
        concept_generator=None, parameter_generator=add_domain_parameters
    )

    sample_2x2_1reward = update_dict(sample_1x3, instances=["sample_2x2_1reward.pddl"])
    sample_2x2_2rewards = update_dict(sample_1x3, instances=["sample_2x2_2rewards.pddl"])
    sample_3x3_2rewards = update_dict(sample_1x3, instances=["sample_3x3_2rewards.pddl"])
    instance_5 = update_dict(sample_1x3, instances=["instance_5.pddl", "instance_4_blocked.pddl"])
    instance_5_no_marking = update_dict(instance_5, complete_only_wrt_optimal=False,)

    parameters = {
        "sample_1x3": sample_1x3,
        "sample_2x2_1reward": sample_2x2_1reward,
        "sample_2x2_2rewards": sample_2x2_2rewards,  # Overfits
        "sample_3x3_2rewards": sample_3x3_2rewards,  # Overfits
        "instance_5": instance_5,
        "instance_5_no_marking": instance_5_no_marking,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def build_expected_features(lang):
    obj_t = lang.Object

    current_cell = PrimitiveConcept(lang.get("at"))
    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    blocked = PrimitiveConcept(lang.get("blocked"))
    unblocked = NotConcept(blocked, obj_t)
    reward = PrimitiveConcept(lang.get("reward"))

    adjacent_unblocked = RestrictRole(adjacent_role, unblocked)

    return [
        ConceptCardinalityFeature(reward),
        # MinDistanceFeature(current_cell, adjacent_role, reward),
        MinDistanceFeature(current_cell, adjacent_unblocked, reward),
    ]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[reward]": "num-rewards",
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
