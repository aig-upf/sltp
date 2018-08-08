#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole, UniversalConcept, GoalConcept, MinDistanceFeature


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "grid-circles"

    domain = "domain.pddl"
    instance = "instance_5.pddl"
    # instance = "sample_1x3.pddl"
    # instance = "sample_2x2_1reward.pddl"
    instance = "sample_2x2_2rewards.pddl"

    steps = generate_pipeline(
        # pipeline="sat",
        pipeline="maxsat",
        domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
        instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

        # Location of the FS planner, used to do the state space sampling
        planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

        # Type of sampling procedure. Only breadth-first search implemented ATM
        driver="bfs",

        # Number of states to be expanded in the sampling procedure
        num_states=604,

        max_concept_size=5,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        # concept_generator=build_expected_concepts,
        # feature_generator=build_expected_features,

        # Whether to use distance features (default: False)
        use_distance_features=True,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=add_domain_parameters,

        feature_namer=feature_namer,

        relax_numeric_increase=True,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # encoding_k=10,
        # encoding_m=10,
    )
    exp = Experiment(steps)
    exp.run()


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


def build_expected_features(lang):
    current_cell = PrimitiveConcept(lang.get("at"))
    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    reward = PrimitiveConcept(lang.get("reward"))
    return [MinDistanceFeature(current_cell, adjacent_role, reward)]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[reward]": "num-uncollected-rewards",
        "bool[And(Not({a}), holding)]": "H",
        "bool[Exists(Inverse(on),{a})]": "Z",
        "card[Exists(Star(on),{a})]": "n(x)",
        "card[And(And(And(Not(Exists(Star(on),{a})), Not(Exists(Star(Inverse(on)),{a}))), Not({a})), Not(holding))]": "m(x)",
    }.get(s, s)


if __name__ == "__main__":
    main()
