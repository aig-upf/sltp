#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import PrimitiveConcept, ConceptCardinalityFeature, MinDistanceFeature, AndConcept, PrimitiveRole, \
    GoalConcept, NominalConcept, ExistsConcept, StarRole, NotConcept


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "taxi"
    domain = "domain.pddl"
    instance = "instance_3.pddl"
    instance = "sample3.pddl"

    experiment1 = dict(instance=instance, num_states=6, max_concept_size=6, max_concept_grammar_iterations=2,
                       distance_feature_max_complexity=5, relax_numeric_increase=True,)

    experiment2 = dict(instance="sample2.pddl", num_states=100, max_concept_size=6, max_concept_grammar_iterations=2,
                       distance_feature_max_complexity=5, relax_numeric_increase=True, )

    defaults = dict(
        pipeline="maxsat",
        # pipeline="sat",

        domain=domain,
        instance=instance,

        # Location of the FS planner, used to do the state space sampling
        planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

        # Type of sampling procedure. Only breadth-first search implemented ATM
        driver="bfs",

        # Number of states to be expanded in the sampling procedure
        num_states=30,

        # Max. size of the generated concepts (mandatory)
        max_concept_size=6,

        # Max. number of iterations of the concept-generation grammar. Optional. Defaults to infinity,
        # in which case the limit is set by max_concept_size alone.
        max_concept_grammar_iterations=2,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        concept_generator=None,

        # Or, alternatively, provide directly the features instead of the concepts (default: None)
        feature_generator=build_expected_features,

        # Max. allowed complexity for distance features (default: 0)
        distance_feature_max_complexity=5,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=None,

        # Use the relaxed (weak) increase semantics
        relax_numeric_increase=False,

        # Optionally, use a method that gives handcrafted names to the features
        # (default: None, which will use their string representation)
        feature_namer=None,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # The number of features and actions for the SAT encoding
        # encoding_k=10,
        # encoding_m=10,
    )

    parameters = {**defaults, **experiment1}  # Copy defaults, overwrite with later dicts
    parameters["domain"] = os.path.join(BENCHMARK_DIR, domain_dir, parameters["domain"])
    parameters["instance"] = os.path.join(BENCHMARK_DIR, domain_dir, parameters["instance"])
    steps = generate_pipeline(**parameters)
    exp = Experiment(steps)
    exp.run()


def add_domain_parameters(language):
    return []


def build_expected_features(lang):
    obj_t = lang.Object
    locp = PrimitiveConcept(lang.get("locp"))
    loct = PrimitiveConcept(lang.get("loct"))
    # loc_fuel = PrimitiveConcept(lang.get("loc_fuel"))
    # min_fuel_level = PrimitiveConcept(lang.get("min_fuel_level"))
    # current_fuel = PrimitiveConcept(lang.get("current_fuel"))
    # destination = GoalConcept(lang.get("locp"))
    destination = GoalConcept(lang.get("locp"))
    inside_taxi = NominalConcept("inside_taxi", obj_t)

    adjacent_role = PrimitiveRole(lang.get("adjacent"))
    # succ = PrimitiveRole(lang.get("succ"))

    # card[Exists(Star(succ),current_fuel)]
    # fuel_amount = ExistsConcept(StarRole(succ), current_fuel)

    # card[Exists(adjacent,And(loct, locp_g))] [k=5, id=36]
    yyy = ExistsConcept(adjacent_role, AndConcept(locp, destination, "cell"))

    return [
        MinDistanceFeature(loct, adjacent_role, locp),  # Distance btw taxi and passenger
        MinDistanceFeature(loct, adjacent_role, destination),  # Distance btw taxi and destination
        ConceptCardinalityFeature(AndConcept(locp, inside_taxi, "cell")),  # Whether passenger is inside the taxi
        ConceptCardinalityFeature(AndConcept(locp, destination, "cell")),  # Whether passenger is at destination

        # MinDistanceFeature(loct, adjacent_role, loc_fuel),  # Distance btw taxi and fuel station
        # ConceptCardinalityFeature(AndConcept(current_fuel, min_fuel_level, "fuel_level")),  # Whether tank is empty
        # ConceptCardinalityFeature(AndConcept(NotConcept(current_fuel, obj_t), min_fuel_level, "fuel_level")),  # Whether tank is full
        # ConceptCardinalityFeature(fuel_amount),
        # ConceptCardinalityFeature(yyy),
    ]


if __name__ == "__main__":
    main()
