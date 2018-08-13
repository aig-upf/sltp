#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "gridworld"

    # domain = "domain.pddl"
    # instance = "instance_3.pddl"
    domain = "domain_strips.pddl"
    instance = "instance_strips_5.pddl"

    parameter_generator = add_domain_parameters_strips if "strips" in domain else add_domain_parameters

    steps = generate_pipeline(pipeline="maxsat",
                              domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

                              # Location of the FS planner, used to do the state space sampling
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

                              # Type of sampling procedure. Only breadth-first search implemented ATM
                              driver="bfs",

                              # Number of states to be expanded in the sampling procedure
                              num_states=60,

                              max_concept_size=10,

                              # Provide a special, handcrafted method to generate concepts, if desired.
                              # This will override the standard concept generation procedure (default: None)
                              # concept_generator=build_expected_concepts,

                              # Max. allowed complexity for distance features (default: 0)
                              distance_feature_max_complexity=10,

                              # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
                              # be used (default: None)
                              parameter_generator=parameter_generator,

                              # What optimization criteria to use in the max-sat problem
                              optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                              # optimization=OptimizationPolicy.NUM_FEATURES
                              )
    exp = Experiment(steps)
    exp.run()


def add_domain_parameters(language):
    # language.constant(2, "coordinate")  # x-goal coordinate
    # language.constant(3, "coordinate")  # x-goal coordinate
    # language.constant(10, "coordinate")  # grid limits!!
    # [language.constant(i, "coordinate") for i in range(1, 11)]
    return [language.constant(1, "coordinate")]  # grid limits!!


def add_domain_parameters_strips(language):
    return []


if __name__ == "__main__":
    main()
