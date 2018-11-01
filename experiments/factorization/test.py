#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment


def experiment(experiment_name=None):
    domain_dir = "blocks"
    domain = "domain.pddl"

    # A small testing instance nonetheless producing an abstraction
    debugging_test = dict(
        instances="instance_4_clear_x.pddl",
        num_states=20, num_sampled_states=10)

    parameters = {
        "test": debugging_test,
    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
