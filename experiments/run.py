#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import importlib
import sys

from defaults import generate_experiment


def report_and_exit(msg):
    print("ERROR: {}".format(msg))
    sys.exit(-1)


if __name__ == "__main__":
    expid = sys.argv[1]

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
    experiment.run(sys.argv[2:])
