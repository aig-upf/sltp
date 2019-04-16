#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import importlib
import sys

from defaults import generate_experiment


def report_and_exit(msg):
    print("ERROR: {}".format(msg))
    sys.exit(-1)


def do(expid, steps):
    name_parts = expid.split(":")
    if len(name_parts) != 2:
        report_and_exit('Wrong experiment ID "{}"'.format(expid))

    scriptname, expname = name_parts

    try:
        mod = importlib.import_module(scriptname)
    except ImportError:
        report_and_exit('No script named "{}".py found on current directory'.format(scriptname))

    try:
        experiments = mod.experiments()
    except AttributeError:
        report_and_exit('Expected method "experiments" not found in script "{}"'.format(expname, scriptname))

    if expname not in experiments:
        report_and_exit('No experiment named "{}" in current experiment script'.format(expname))

    parameters = experiments[expname]
    experiment = generate_experiment(**parameters)
    experiment.run(steps)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE: {} domain:experiment_name <steps>,\n\te.g.:\t./run.py logistics:p1 --all".format(sys.argv[0]))
        sys.exit(-1)
    do(sys.argv[1], sys.argv[2:])
