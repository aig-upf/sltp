#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from abstractions_defaults import generate_experiment


def experiment(experiment_name=None):
    domain_dir = "logistics98"
    domain = "domain.pddl"

    sample1 = dict(
        instance="sample1.pddl",
        num_states=100, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    sample1_rnd = dict(
        instance="sample1.pddl",
        num_states=1000, num_sampled_states=60, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    sample2 = sample1.copy()
    sample2.update(dict(instance="sample2.pddl"))

    sample2_rnd = sample1_rnd.copy()
    sample2_rnd.update(dict(instance="sample2.pddl"))

    parameters = {
        "sample1": sample1,
        "sample1_rnd": sample1_rnd,
        "sample2": sample2,
        "sample2_rnd": sample2_rnd,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
