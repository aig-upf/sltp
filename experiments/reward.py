#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from sltp.driver import run_experiment

from sltp.util.defaults import generate_experiment
from sltp.util.misc import update_dict
from tarski.dl import AndConcept, PrimitiveConcept, NotConcept, PrimitiveRole, UniversalConcept, MinDistanceFeature, ConceptCardinalityFeature, RestrictRole


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

    exps = dict()

    # One reason for overfitting: in a 3x3 grid, with 2 booleans per dimension you can perfectly represent any position
    exps['sample_1x3'] = dict(
        instances=["sample_1x3.pddl"],
        complete_only_wrt_optimal=True,
        num_states=1500, max_concept_size=8, max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        feature_namer=feature_namer,
        # feature_generator=build_expected_features,
        concept_generator=None, parameter_generator=add_domain_parameters
    )

    exps['sample_2x2_1reward'] = update_dict(exps['sample_1x3'], instances=["sample_2x2_1reward.pddl"])
    exps['sample_2x2_2rewards'] = update_dict(exps['sample_1x3'], instances=["sample_2x2_2rewards.pddl"])
    exps['sample_3x3_2rewards'] = update_dict(exps['sample_1x3'], instances=["sample_3x3_2rewards.pddl"])
    exps['instance_5'] = update_dict(exps['sample_1x3'], instances=["instance_5.pddl", "instance_4_blocked.pddl"])
    exps['instance_5_no_marking'] = update_dict(exps['instance_5'], complete_only_wrt_optimal=False,)

    # Same but using goal-concepts instead of goal parameters:
    exps["instance_5_gc"] = update_dict(exps["instance_5"], parameter_generator=None)

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
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
    run_experiment(exp, sys.argv[2:])

