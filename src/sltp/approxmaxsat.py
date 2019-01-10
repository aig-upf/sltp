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
import logging
import re

from sltp.errors import CriticalPipelineError
from sltp.util.command import execute

from .solvers import solve
from .returncodes import ExitCode


def run(config, data, rng):
    logging.info("Translating feature matrix into WSAT format".format())

    # We assume that blai's translator and solver are on the path and called 'wsat-translator' and  'wsat-solver'
    command = "wsat-translator {}".format(config.sat_theory_prefix)
    # retcode = execute(command=command.split(' '),
                      # stdout=config.sample_file,
                      # cwd=config.planner_location
                      # )
    # if retcode:
    #     raise CriticalPipelineError("Error running MAXSAT translator. Check the logs. Return code: {}".format(retcode))

    command = "wsat-solver --implicit {}".format(config.sat_theory_prefix)
    # command = "wsat-solver /home/frances/projects/code/weighted-max-sat/examples/blocks04_theory.wsat".format()

    retcode = execute(command=command.split(' '),
                      stdout=config.maxsat_solver_out,
                      )
    if retcode:
        raise CriticalPipelineError("Error running MAXSAT solver. Check the logs. Return code: {}".format(retcode))

    # Parse the output
    with open(config.maxsat_solver_out, 'r') as f:
        out = f.read()

    # We parse an output like "solution=[0,16,17,2,21]"
    try:
        start = out.find('solution=[') + 10
        end = out.find(']', start)
        features_str = out[start:end]
        # features_str = re.search('solution=[(.+?)]', out).group(1)
    except AttributeError:
        raise CriticalPipelineError("No solution found in output of MAXSAT solver, in file '{}'".format(
            config.maxsat_solver_out))

    features = list(map(int, features_str.split(',')))

    # TODO Check output when the problem is UNSAT
    # if unsat:
    #     return ExitCode.MaxsatModelUnsat, dict()

    return ExitCode.Success, dict(selected_feature_idxs=features)


