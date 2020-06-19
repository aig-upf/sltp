#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.util.defaults import generate_experiment


def experiment(experiment_name=None):
    domain_dir = "gripper"
    domain = "domain.pddl"

    exps = dict()

    exps["sample_small"] = dict(
        instances="sample-small.pddl",
        test_domain=domain,
        test_instances=["prob03.pddl", "prob04.pddl"],
        #test_instances=[],
        num_states=500,
        max_concept_size=10,
        max_concept_grammar_iterations=3,
        concept_generator=None,
        parameter_generator=add_domain_parameters,
        #pipeline="maxsat_poly",
        feature_namer=feature_namer,)

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    return generate_experiment(domain_dir, domain, **exps[experiment_name])


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
        "card[Exists(at,{roomb})]": "nballs-B",
        "card[Exists(carry,<universe>)]": "ncarried",
        "bool[And(at-robby, {roomb})]": "robot-at-B",
        "card[Exists(at,Not(at-robby))]": "nballs-in-rooms-with-no-robot",
        "card[free]": "nfree-grippers",
        "bool[Exists(at-robby,{roomb})]": "robby-is-at-B",
        "card[Exists(at,Exists(Inverse(at-robby),<universe>))]": "nballs-in-room-with-some-robot",
        "card[And(Exists(gripper,Exists(at-robby,{roomb})),free)]": "nfree-grippers-at-B",
        "card[Exists(at-robby,{roomb})]": "nrobots-at-B",
        "card[Exists(gripper,Exists(at-robby,{roomb}))]": "ngrippers-at-B",
        "card[Exists(carry,Exists(gripper,Exists(at-robby,{roomb})))]": "nballs-carried-in-B",
        "card[Exists(at,And(Forall(Inverse(at-robby),<empty>), Not({roomb})))]": "nballs-in-some-room-notB-without-any-robot",
        "bool[And(Exists(Inverse(at),<universe>), And({roomb}, Not(at-robby)))]": "some-ball-in-B-but-robbot-not-in-B",
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
