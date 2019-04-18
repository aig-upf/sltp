#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import importlib
import os
import sys

from sltp.util import console
from sltp.util.bootstrap import setup_argparser
from sltp.util.misc import update_dict

from sltp.util.defaults import generate_experiment


def report_and_exit(msg):
    print("ERROR: {}".format(msg))
    sys.exit(-1)


def experiment(experiment_name=None, **kwargs):

    base = update_dict(kwargs,
                       domain_dir="gripper",
                       domain="domain.pddl")

    exps = dict()
    exps["aaai_prob01"] = update_dict(
        base,
        instances=["prob01.pddl", "sample02.pddl"],
        num_states=2000, max_width=[-1],
        num_sampled_states=100,
        complete_only_wrt_optimal=True,
        max_concept_size=8, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    return generate_experiment(**exps[experiment_name])


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
        "card[Exists(at,{roomb})]": "nballs-B",
    }.get(s, s)


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def main():
    args = setup_argparser().parse_args(sys.argv[1:])

    if args.workspace is not None:
        exp = experiment(args.exp_id, workspace=os.path.abspath(args.workspace))
    else:
        exp = experiment(args.exp_id)

    print(args)

    exp.run(args.steps)


if __name__ == "__main__":
    main()
