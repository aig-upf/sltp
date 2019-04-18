#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from tarski.dl import PrimitiveConcept, GoalConcept, PrimitiveRole

from sltp.util.defaults import generate_experiment


def experiment(experiment_name=None):
    domain_dir = "taxi"
    domain = "domain.pddl"

    exps = dict()

    exps["simple"] = dict(
        instances="instance_3.pddl",
        test_domain=domain, test_instances=["instance_5.pddl"],
        num_states=300,
        max_concept_size=10,
        distance_feature_max_complexity=5,
        parameter_generator=None,
        feature_namer=feature_namer,
        concept_generator=build_expected_concepts,
        complete_only_wrt_optimal=True
    )

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    return generate_experiment(domain_dir, domain, **exps[experiment_name])


def feature_namer(feature):
    s = str(feature)
    return {
    }.get(s, s)


def build_expected_concepts(lang):
    """  """
    # obj_t = lang.Object
    loct, locp, adj = lang.get("loct", "locp", "adjacent")
    concepts = [PrimitiveConcept(loct), PrimitiveConcept(locp), GoalConcept(locp)]
    roles = [PrimitiveRole(adj)]
    return [], concepts, roles  # atoms, concepts, roles


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
