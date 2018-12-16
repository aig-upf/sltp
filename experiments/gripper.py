#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment
from sltp.util.misc import update_dict


def experiment(experiment_name=None):
    domain_dir = "gripper"
    # domain_dir = "gripper-m"
    domain = "domain.pddl"

    sample_small = dict(
        instance="sample-small.pddl",
        num_states=200, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    prob01 = dict(
        instance="prob01.pddl",
        num_states=300, num_sampled_states=None, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    #
    prob01_rnd = dict(
        instance="prob01.pddl",
        num_states=2000, num_sampled_states=50, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    #
    aaai_prob01 = dict(
        instances=["prob01.pddl", "prob02.pddl"],
        num_states=2000, max_width=[-1],
        num_sampled_states=100,
        complete_only_wrt_optimal=True,
        max_concept_size=8, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    aaai_prob01_no_marking = update_dict(aaai_prob01, complete_only_wrt_optimal=False)

    parameters = {
        "sample_small": sample_small,
        "prob01": prob01,
        "prob01_rnd": prob01_rnd,
        "aaai_prob01": aaai_prob01,
        "aaai_prob01_no_marking": aaai_prob01_no_marking

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
        "card[Exists(at,{roomb})]": "nballs-B",
        "card[Exists(carry,<universe>)]": "ncarried",
        "bool[And(at-robby, {roomb})]": "robot-at-B",
        "card[free]": "nfree-grippers",
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
