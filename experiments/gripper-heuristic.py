#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from defaults import generate_experiment
from tarski.dl import PrimitiveRole, NominalConcept, ExistsConcept, NotConcept, UniversalConcept, AndConcept, \
    ForallConcept, EmptyConcept


def experiment(experiment_name=None):
    domain_dir = "gripper-m"
    domain = "domain.pddl"

    prob01 = dict(
        pipeline="heuristic",
        lp_max_weight=10,
        instances="prob01.pddl",
        num_states=300, num_sampled_states=None, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    parameters = {
        "prob01": prob01,
    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def generate_chosen_concepts(lang):
    """  """
    # card[Exists(at,Not({roomb}))] 4    C1
    # card[Exists(carry,<universe>)] 2   C2
    # bool[Exists(at-robby,{roomb})] 3   RX
    # card[Exists(gripper,Exists(at-robby,{roomb}))] 5   (intermediate)
    # card[Exists(carry,Exists(gripper,Exists(at-robby,{roomb})))] 7

    obj_t = lang.Object

    at = PrimitiveRole(lang.get("at"))
    carry = PrimitiveRole(lang.get("carry"))
    at_robby = PrimitiveRole(lang.get("at-robby"))
    gripper = PrimitiveRole(lang.get("gripper"))
    x_param = NominalConcept("roomb", obj_t)
    c1 = ExistsConcept(at, NotConcept(x_param, obj_t))
    c2 = ExistsConcept(carry, UniversalConcept("object"))
    rx = ExistsConcept(at_robby, x_param)
    c3 = ExistsConcept(gripper, rx)
    c4 = ExistsConcept(carry, c3)

    concepts = [c1, c2, c3, c4]
    return [], concepts, []  # atoms, concepts, roles


def debug_weird_concept(lang):
    #  card[And(Forall(carry,<empty>), Forall(at-robby,{roomb}))]

    obj_t = lang.Object
    bot = EmptyConcept("object")

    carry = PrimitiveRole(lang.get("carry"))
    at_robby = PrimitiveRole(lang.get("at-robby"))
    x_param = NominalConcept("roomb", obj_t)

    c1 = ForallConcept(carry, bot)
    c2 = ForallConcept(at_robby, x_param)
    c = AndConcept(c1, c2, "object")

    concepts = [c]
    return [], concepts, []  # atoms, concepts, roles


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
    }.get(s, s)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
