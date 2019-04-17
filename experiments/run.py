#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import importlib
import sys

from sltp.util.bootstrap import setup_argparser

from defaults import generate_experiment


def report_and_exit(msg):
    print("ERROR: {}".format(msg))
    sys.exit(-1)


def do(expid, steps=[], workspace=None):
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

    # Recreate the namespace, that will likely be parsed again
    args_ = [expid]
    args_ = args_ if not steps else args_ + [str(s) for s in steps]
    args_ = args_ if workspace is None else args_ + ["--workspace", workspace]
    experiment.run(args_)


if __name__ == "__main__":
    args = setup_argparser().parse_args(sys.argv[1:])
    do(args.__dict__["domain:experiment"], args.steps, args.workspace)
