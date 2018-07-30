#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_full_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "gridworld"

    # domain = "domain.pddl"
    # instance = "instance_3.pddl"
    domain = "domain_strips.pddl"
    instance = "instance_strips_5.pddl"

    steps = generate_full_pipeline(domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
                                   instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),
                                   driver="bfs",
                                   planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                                   num_states=60,
                                   concept_depth=2,
                                   optimization=OptimizationPolicy.TOTAL_FEATURE_DEPTH,
                                   use_distance_features=True
                                   # optimization=OptimizationPolicy.NUM_FEATURES
                                   )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
