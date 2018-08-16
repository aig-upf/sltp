#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from abstractions_defaults import generate_experiment
from common import build_ijcai_paper_bw_concepts, ijcai_paper_bw_feature_namer, no_parameter


def experiment(experiment_name=None):
    domain_dir = "blocks-tower"
    domain = "domain.pddl"

    sample_4 = dict(
        instance="sample_4.pddl",
        num_states=200, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=no_parameter,
        feature_namer=ijcai_paper_bw_feature_namer,)

    #
    sample_5_rnd = dict(
        instance="sample_5.pddl",
        num_states=2000, num_sampled_states=40, random_seed=12,
        max_concept_size=15, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=no_parameter,
        feature_namer=ijcai_paper_bw_feature_namer,)

    parameters = {
        "sample_4": sample_4,
        "sample_5": sample_5_rnd,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
