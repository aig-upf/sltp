#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from experiments.common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, ijcai_paper_bw_feature_namer


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "blocks"
    domain = "domain.pddl"
    # instance = "probBLOCKS-4-0.pddl"
    # instance = "instance_3_clear_x.pddl"
    instance = "instance_4_clear_x.pddl"
    # instance = "instance_5_clear_x.pddl"

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
        num_states=80,

        # Max. size of the generated concepts (mandatory)
        max_concept_size=8,

        # Max. number of iterations of the concept-generation grammar. Optional. Defaults to infinity,
        # in which case the limit is set by max_concept_size alone.
        max_concept_grammar_iterations=3,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        # concept_generator=build_ijcai_paper_bw_concepts,

        # Or, alternatively, provide directly the features instead of the concepts (default: None)
        feature_generator=None,

        # Whether to use distance features (default: False)
        # use_distance_features=True,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=add_bw_domain_parameters,

        # Use the relaxed (weak) increase semantics
        relax_numeric_increase=False,

        # Optionally, use a method that gives handcrafted names to the features
        # (default: None, which will use their string representation)
        feature_namer=ijcai_paper_bw_feature_namer,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # The number of features and actions for the SAT encoding
        encoding_k=10,
        encoding_m=10,
    )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
