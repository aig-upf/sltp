#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from abstractions_defaults import generate_experiment
from common import update_dict


def experiment(experiment_name=None):
    domain_dir = "gridworld"
    # domain_dir = "gripper-m"
    domain = "domain_strips.pddl"

    # This one overfits: in a 3x3 grid, with 2 booleans per dimension you can represent
    # the position
    sample_3 = dict(
        instance="instance_strips_3.pddl",
        num_states=80, max_concept_size=10, max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        concept_generator=None, parameter_generator=add_domain_parameters_strips)

    sample_5 = update_dict(sample_3, instance="instance_strips_5.pddl", num_states=100)
    sample_5_not_corner = update_dict(sample_5, instance="instance_strips_5_not_corner.pddl")
    sample_5_rnd = update_dict(sample_5, num_states=1000, num_sampled_states=80, random_seed=12)
    sample_5_not_corner_rnd = update_dict(sample_5_not_corner, num_states=1000, num_sampled_states=80, random_seed=12)

    parameters = {
        "sample_3": sample_3,
        "sample_5": sample_5,
        "sample_5_not_corner": sample_5_not_corner,
        "sample_5_rnd": sample_5_rnd,
        "sample_5_not_corner_rnd": sample_5_not_corner_rnd,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def add_domain_parameters(language):
    # language.constant(2, "coordinate")  # x-goal coordinate
    # language.constant(3, "coordinate")  # x-goal coordinate
    # language.constant(10, "coordinate")  # grid limits!!
    # [language.constant(i, "coordinate") for i in range(1, 11)]
    return [language.constant(1, "coordinate")]  # grid limits!!


def add_domain_parameters_strips(language):
    return []


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
