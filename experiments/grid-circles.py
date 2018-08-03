#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "grid-circles"

    domain = "domain_strips.pddl"
    instance = "instance_strips_3.pddl"

    steps = generate_pipeline(pipeline="maxsat",
                              domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

                              # Location of the FS planner, used to do the state space sampling
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

                              # Type of sampling procedure. Only breadth-first search implemented ATM
                              driver="bfs",

                              # Number of states to be expanded in the sampling procedure
                              num_states=1000,

                              # Number of iterations of the concept generation grammar.
                              concept_depth=2,

                              # Provide a special, handcrafted method to generate concepts, if desired.
                              # This will override the standard concept generation procedure (default: None)
                              concept_generator=build_expected_concepts,

                              # Whether to use distance features (default: False)
                              use_distance_features=True,

                              # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
                              # be used (default: None)
                              parameter_generator=add_domain_parameters,

                              # What optimization criteria to use in the max-sat problem
                              optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                              # optimization=OptimizationPolicy.NUM_FEATURES
                              )
    exp = Experiment(steps)
    exp.run()


def add_domain_parameters(language):
    return []


def build_expected_concepts(lang):
    obj_t = lang.get_sort("cell")

    current_cell = PrimitiveConcept(lang.get("at"))
    reward_cells = PrimitiveConcept(lang.get("reward"))
    blocked_cells = PrimitiveConcept(lang.get("blocked"))
    unblocked_cells = NotConcept(blocked_cells, obj_t)
    visited_cells = PrimitiveConcept(lang.get("visited"))
    unvisited_cells = NotConcept(visited_cells, obj_t)
    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    unvisited_cells_with_reward = AndConcept(unvisited_cells, reward_cells, "cell")

    concepts = [current_cell, unvisited_cells, unvisited_cells_with_reward, unblocked_cells]
    roles = [adjacent_role]

    return [], concepts, roles  # atoms, concepts, roles


if __name__ == "__main__":
    main()
