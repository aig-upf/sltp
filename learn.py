#!/usr/bin/env python

#  Copyright (C) 2018-<date> Blai Bonet
#
#  Permission is hereby granted to distribute this software for
#  non-commercial research purposes, provided that this copyright
#  notice is included with any such distribution.
#
#  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
#  EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
#  SOFTWARE IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
#  ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#  Blai Bonet, bonet@ldc.usb.ve, bonetblai@gmail.com
#  Guillem Frances, guillem.frances@unibas.ch
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

from util.bootstrap import bootstrap
from util.console import print_header
import features as fgenerator
from solvers import solve

signal(SIGPIPE, SIG_DFL)

BASEDIR = os.path.dirname(os.path.realpath(__file__))
VERSION = "0.1"


def main(args):

    all_features, states, goal_states, transitions, cache = fgenerator.main(args)

    rundir = os.path.join(BASEDIR, "runs")
    cnf_filename = os.path.join(rundir, args.result_filename + ".cnf")

    translator = ModelTranslator(all_features, states, goal_states, transitions, cache)
    translator.run(cnf_filename)

    # solution = solve(rundir, cnf_filename, 'wpm3')
    # solution = solve(rundir, cnf_filename, 'maxino')
    solution = solve(rundir, cnf_filename, 'openwbo')

    if not solution.solved and solution.result == "UNSATISFIABLE":
        print("The MAXSAT encoding is UNSATISFIABLE")
    else:
        translator.decode_solution(solution.assignment)


def run_fs_planner():
    pass
'./run.py --instance /home/frances/projects/code/downward-benchmarks/blocks/probBLOCKS-4-0.pddl  --driver bfs --disable-static-analysis --options="max_expansions=100" > blocks.txt'


def run_feature_generator():
    pass


def run_maxsat_encoder():
    pass


def run_maxsat_solver():
    pass


def run_action_model_generator():
    pass


if __name__ == "__main__":
    print_header("Generalized Action Model Learner, v.{}".format(VERSION))
    _args = bootstrap(sys.argv[1:])
    main(_args)


