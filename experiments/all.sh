#!/bin/bash

set -e -x

#STEPS=""  # Run all pipeline steps
STEPS="1 2 3"  # Generate CNF input matrices only

./run.py blocks:aaai_clear_x_simple_hybrid ${STEPS}
./run.py blocks:aaai_bw_on_x_y_completeness_opt ${STEPS}
./run.py gripper:aaai_prob01 ${STEPS}
./run.py logistics:p1 ${STEPS}
./run.py barman:p1 ${STEPS}
./run.py childsnack:p1 ${STEPS}
./run.py spanner:p1 ${STEPS}
./run.py grid:p1 ${STEPS}
./run.py miconic:p1 ${STEPS}
./run.py hanoi:p1 ${STEPS}
./run.py visitall:p1 ${STEPS}
./run.py depot:p1 ${STEPS}
./run.py satellite:p1 ${STEPS}
./run.py blocks2:arbitrary1 ${STEPS}
./run.py blocks2:on_x_y ${STEPS}
./run.py blocks2:on_x_y_gc ${STEPS}