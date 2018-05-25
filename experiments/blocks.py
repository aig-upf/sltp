#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os

from learn import Experiment, generate_full_pipeline, BENCHMARK_DIR


def main():

    steps = generate_full_pipeline(domain=os.path.join(BENCHMARK_DIR, "blocks", "domain.pddl"),
                                   instance=os.path.join(BENCHMARK_DIR, "blocks", "instance_3_clear_x.pddl"),
                                   driver="bfs",
                                   planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                                   num_states=30,
                                   concept_depth=2)
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
