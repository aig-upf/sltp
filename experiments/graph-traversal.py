#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole, UniversalConcept, GoalConcept, MinDistanceFeature, ConceptCardinalityFeature, \
    EmpiricalBinaryConcept


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "graph-traversal"

    domain = "domain.pddl"
    instance = "instance_5_12.pddl"
    instance = "instance_10_50.pddl"

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
        num_states=600,

        max_concept_size=5,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        # concept_generator=build_expected_concepts,
        # feature_generator=build_expected_features,

        # Whether to use distance features (default: False)
        use_distance_features=True,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        # parameter_generator=add_domain_parameters,

        feature_namer=feature_namer,

        relax_numeric_increase=True,
        # relax_numeric_increase=False,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES
    )
    exp = Experiment(steps)
    exp.run()


def add_domain_parameters(language):
    return []


def build_expected_concepts(lang):
    return [], [], []  # atoms, concepts, roles


def build_expected_features(lang):
    current_cell = PrimitiveConcept(lang.get("at"))
    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    at_g = GoalConcept(lang.get("at"))
    return [
        MinDistanceFeature(current_cell, adjacent_role, at_g),
        ConceptCardinalityFeature(AndConcept(current_cell, at_g, "cell"))
    ]


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


if __name__ == "__main__":
    main()
