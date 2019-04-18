#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import importlib
import sys

from sltp.util import console
from sltp.util.bootstrap import setup_argparser

from defaults import generate_experiment


def report_and_exit(msg):
    print("ERROR: {}".format(msg))
    sys.exit(-1)


def do(expid, steps=None, workspace=None, show_steps_only=False):
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
    if workspace is not None:
        parameters["workspace"] = workspace

    experiment = generate_experiment(**parameters)

    if show_steps_only:
        console.print_hello()
        print('Experiment with id "{}" is configured with the following steps:'.format(expid))
        print(experiment.print_description())
        return

    experiment.run(steps)


def main():
    args = setup_argparser().parse_args(sys.argv[1:])
    do(args.exp_id, args.steps, args.workspace, args.show)


if __name__ == "__main__":
    main()
