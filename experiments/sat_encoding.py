#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_full_pipeline, BENCHMARK_DIR, generate_alt_pipeline

    domain_dir = "blocks"
    steps = generate_alt_pipeline(domain=os.path.join(BENCHMARK_DIR, domain_dir, "domain.pddl"),
                                  instance=os.path.join(BENCHMARK_DIR, domain_dir, "instance_5_clear_x.pddl"),
                                  driver="bfs",
                                  planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                                  num_states=30,
                                  concept_depth=1,
                                  encoding_k=10,
                                  encoding_m=10,
                                  use_distance_features=False
                                  )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
