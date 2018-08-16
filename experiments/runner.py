#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

import blocks
import gripper

DOMAINS = {
    "blocks": blocks,
    "gripper": gripper,
}


def create_parser():
    parser = argparse.ArgumentParser(description='Let\'s run some experiments!')
    parser.add_argument('-d', '--domain', required=True, help="Which domain to run the experiment on.", choices=DOMAINS.keys())
    parser.add_argument('-e', '--experiment', required=True, help="Experiment name.")
    return parser


def run(args):
    parser = create_parser()
    args = parser.parse_args(args)
    domain = DOMAINS[args.domain]
    exp = domain.experiment(args.experiment)
    exp.run(["--all"])


if __name__ == "__main__":
    run(sys.argv[1:])