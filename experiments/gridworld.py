#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from abstractions_defaults import generate_experiment


def experiment(experiment_name=None):
    domain_dir = "gridworld"
    # domain_dir = "gripper-m"
    domain = "domain_strips.pddl"

    sample_3 = dict(
        instance="instance_strips_3.pddl",
        num_states=80, max_concept_size=10, max_concept_grammar_iterations=3,
        distance_feature_max_complexity=8,
        concept_generator=None, parameter_generator=add_domain_parameters_strips)

    sample_5 = sample_3.copy()
    sample_5.update(dict(instance="instance_strips_5.pddl"), num_states=100)

    sample_5_rnd = sample_5.copy()
    sample_5_rnd.update(dict( num_states=1000, num_sampled_states=80, random_seed=12))

    parameters = {
        "sample_3": sample_3,
        "sample_5": sample_5,
        "sample_5_rnd": sample_5_rnd,

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
