#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from tarski.dl import ForallConcept, StarRole, NotConcept, PrimitiveRole, NominalConcept, InverseRole, AndConcept

from abstractions_defaults import generate_experiment
from common import build_ijcai_paper_bw_concepts, ijcai_paper_bw_feature_namer, no_parameter


def experiment(experiment_name=None):
    domain_dir = "blocks-tower"
    domain = "domain.pddl"

    sample_4_1 = dict(
        instance="sample_4.pddl",
        num_states=200, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=no_parameter,
        feature_namer=ijcai_paper_bw_feature_namer,)

    sample_4_2 = sample_4_1.copy()
    sample_4_2.update(dict(max_concept_size=15))

    sample_5_base = dict(
        instance="sample_5_base_x.pddl",
        num_states=2000,
        num_sampled_states=50, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        # concept_generator=build_x_tower_concepts,
        parameter_generator=one_block_fs,
        feature_namer=ijcai_paper_bw_feature_namer,)

    sample_5_base_2 = sample_5_base.copy()
    sample_5_base_2.update(dict(num_sampled_states=60, random_seed=10))

    # This one overfits!
    sample_5_rnd = dict(
        instance="sample_5.pddl",
        num_states=2000, num_sampled_states=40, random_seed=12,
        max_concept_size=15, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=no_parameter,
        feature_namer=ijcai_paper_bw_feature_namer,)

    ssample_5_rnd_2 = sample_5_rnd.copy()
    ssample_5_rnd_2.update(dict(num_sampled_states=70, random_seed=812))

    parameters = {
        "sample_4_1": sample_4_1,
        "sample_4_2": sample_4_2,
        "sample_5_rnd": sample_5_rnd,
        "ssample_5_rnd_2": ssample_5_rnd_2,
        "sample_5_base": sample_5_base,
        "sample_5_base_2": sample_5_base_2,

    }.get(experiment_name or "test")

    return generate_experiment(domain_dir, domain, **parameters)


def one_block_fs(language):
    return [language.constant("b1", "block")]


def build_x_tower_concepts(lang):
    obj_t = lang.Object
    x_nominal = NominalConcept("b1", obj_t)
    table_nominal = NominalConcept("table", obj_t)

    not_x = NotConcept(x_nominal, obj_t)
    not_table = NotConcept(table_nominal, obj_t)

    on_r = PrimitiveRole(lang.get("loc"))
    above = StarRole(on_r)
    below = StarRole(InverseRole(on_r))
    all_below_not_x = ForallConcept(above, not_x)

    not_x_and_all_below_not_x = AndConcept(not_table, AndConcept(all_below_not_x, not_x, "object"), "object")

    concepts = [not_x_and_all_below_not_x]
    return [], concepts, []  # atoms, concepts, roles




if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
