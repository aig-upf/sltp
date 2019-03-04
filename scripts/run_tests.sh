#!/bin/bash

set -e -x

# PYTHONPATH=$PYTHONPATH:src pytest
cd experiments
./gripper.py aaai_prob01_gc --all
./blocks.py aaai_clear_x_simple_hybrid --all
./blocks.py aaai_bw_on_x_y_completeness_opt --all
#./reward.py instance_5 --all

