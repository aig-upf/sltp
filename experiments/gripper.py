#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "gripper"
    # domain_dir = "gripper-m"
    domain = "domain.pddl"
    instance = "prob01.pddl"

    steps = generate_pipeline(
        pipeline="maxsat",
        # pipeline="sat",
        domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
        instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

        # Location of the FS planner, used to do the state space sampling
        planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

        # Type of sampling procedure. Only breadth-first search implemented ATM
        driver="bfs",

        # Number of states to be expanded in the sampling procedure
        num_states=250,

        # Max. size of the generated concepts (mandatory)
        max_concept_size=10,

        # Max. number of iterations of the concept-generation grammar. Optional. Defaults to infinity,
        # in which case the limit is set by max_concept_size alone.
        max_concept_grammar_iterations=3,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        # concept_generator=build_ijcai_paper_bw_concepts,

        # Or, alternatively, provide directly the features instead of the concepts (default: None)
        feature_generator=None,

        # Max. allowed complexity for distance features (default: 0)
        # distance_feature_max_complexity=10,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=add_domain_parameters,

        # Use the relaxed (weak) increase semantics
        relax_numeric_increase=False,

        # Optionally, use a method that gives handcrafted names to the features
        # (default: None, which will use their string representation)
        feature_namer=feature_namer,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # The number of features and actions for the SAT encoding
        # encoding_k=10,
        # encoding_m=10,
    )
    exp = Experiment(steps)
    exp.run()


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
        "card[Exists(at,{roomb})]": "nballs-B",
        "card[Exists(carry,<universe>)]": "ncarried",
        "bool[And(at-robby, {roomb})]": "robot-at-B",
        "card[free]": "nfree-grippers",
    }.get(s, s)


if __name__ == "__main__":
    main()
