#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "gripper"
    steps = generate_pipeline(domain=os.path.join(BENCHMARK_DIR, domain_dir, "domain.pddl"),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, "prob01.pddl"),
                              driver="bfs",
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                              num_states=250,
                              concept_depth=1,
                              optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY
                              # optimization=OptimizationPolicy.NUM_FEATURES
                              )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
