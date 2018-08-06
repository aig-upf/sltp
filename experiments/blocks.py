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
    instance = "instance_5_clear_x.pddl"

    steps = generate_pipeline(pipeline="maxsat",
                              domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

                              # Location of the FS planner, used to do the state space sampling
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

                              # Type of sampling procedure. Only breadth-first search implemented ATM
                              driver="bfs",

                              # Number of states to be expanded in the sampling procedure
                              num_states=90,

                              # Maximum depth of the generated concepts
                              max_concept_size=10,

                              # Provide a special, handcrafted method to generate concepts, if desired.
                              # This will override the standard concept generation procedure (default: None)
                              concept_generator=build_ijcai_paper_bw_concepts,

                              # Whether to use distance features (default: False)
                              # use_distance_features=True,

                              # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
                              # be used (default: None)
                              parameter_generator=add_bw_domain_parameters,

                              # Optionally, use a method that gives handcrafted names to the features
                              # (default: use their string representation)
                              feature_namer=ijcai_paper_bw_feature_namer,

                              # What optimization criteria to use in the max-sat problem
                              optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                              # optimization=OptimizationPolicy.NUM_FEATURES
                              )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
