#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from experiments.common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR

    domain_dir = "blocks"
    steps = generate_pipeline(pipeline="sat",
                              domain=os.path.join(BENCHMARK_DIR, domain_dir, "domain.pddl"),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, "instance_5_clear_x.pddl"),
                              driver="bfs",
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                              num_states=90,
                              max_concept_size=10,
                              encoding_k=5,
                              encoding_m=8,
                              concept_generator=build_ijcai_paper_bw_concepts,
                              parameter_generator=add_bw_domain_parameters,
                              # Max. allowed complexity for distance features (default: 0)
                              # distance_feature_max_complexity=10,
                              )
    exp = Experiment(steps)
    exp.run()


if __name__ == "__main__":
    main()
