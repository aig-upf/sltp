#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from experiments.abstractions_defaults import generate_experiment
from experiments.common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, ijcai_paper_bw_feature_namer, \
    build_alt_feature_set
from util.serialization import deserialize_from_string


def main():
    import sys
    sys.path.insert(0, '..')

    domain_dir = "blocks"
    domain = "domain.pddl"
    # instance = "probBLOCKS-4-0.pddl"
    # instance = "instance_4_clear_x.pddl"
    # instance = "instance_5_clear_x.pddl"

    # concept_generator = build_ijcai_paper_bw_concepts
    concept_generator = None

    # A small testing instance nonetheless producing an abstraction
    debugging_test = dict(
        instance="instance_4_clear_x.pddl",
        num_states=20, num_sampled_states=10,
        max_concept_size=3, max_concept_grammar_iterations=1,
        parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # Learns a simple action model which is however overfit to 3 blocks,
    # and not sound in general
    simple_clear_3 = dict(
        instance="instance_3_clear_x.pddl",
        num_states=100, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=concept_generator, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # This example shows that with 4 blocks, even if we expand all states, the model is still overfit to the 4 blocks
    simple_clear_4 = dict(
        instance="instance_4_clear_x.pddl",
        num_states=150, max_concept_size=10, max_concept_grammar_iterations=3, num_sampled_states=None,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # With these settings we generate the desired m(x):
    # card[And(And(Forall(Star(on),Not({a})), Forall(Star(Inverse(on)),Not({a}))), And(Not(holding), Not({a})))] 18
    # And indeed learnt the correct state space!!!
    we_learn_ijcai_features_on_clear_5 = dict(
        instance="instance_5_clear_x_1.pddl",
        num_states=2000, num_sampled_states=40, random_seed=12,
        max_concept_size=18, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    check_ijcai_features_on_clear_5 = dict(
        instance="instance_5_clear_x.pddl",
        num_states=1000, max_concept_size=1, max_concept_grammar_iterations=1, num_sampled_states=100, random_seed=2,
        concept_generator=build_ijcai_paper_bw_concepts, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    check_clear4_features_on_clear_5 = dict(
        instance="instance_5_clear_x.pddl",
        num_states=600, max_concept_size=1, max_concept_grammar_iterations=1,
        concept_generator=build_ijcai_paper_bw_concepts, parameter_generator=add_bw_domain_parameters,
        feature_generator=try_clear5_features,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exp = generate_experiment(domain_dir, domain, **we_learn_ijcai_features_on_clear_5)
    # exp = generate_experiment(domain_dir, domain, **check_clear4_features_on_clear_5)
    exp.run()


def try_clear4_features(lang):
    """ A model learnt by generating the whole 4-block state space """
    features = '[{"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -1277705213948677935, "name": "holding", "size": 1, "sort": "object"}, "hash": -2841485854126214806}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.UniversalConcept", "hash": -9223372036852258203, "size": 0, "sort": "object"}, "hash": -7885522037927662577, "r": {"py/object": "tarski.dl.concepts.PrimitiveRole", "hash": -3816393078228064225, "name": "on", "size": 1, "sort": ["object", "object"]}, "size": 2, "sort": "object"}, "hash": 911132480848590796}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NominalConcept", "hash": -7401258638559185829, "name": "a", "size": 1, "sort": "object"}, "hash": 1733562239334784851, "r": {"py/object": "tarski.dl.concepts.StarRole", "hash": 2196893821251177549, "r": {"py/id": 6}, "size": 2, "sort": {"py/id": 7}}, "size": 4, "sort": "object"}, "hash": 8416726960899686352}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NotConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 5962884402008325778, "name": "ontable", "size": 1, "sort": "object"}, "hash": 2613871271492961489, "size": 2, "sort": "object"}, "hash": 3123638717406933755, "r": {"py/id": 6}, "size": 4, "sort": "object"}, "hash": -8435674864905251064}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 3655234103021149859, "name": "clear", "size": 1, "sort": "object"}, "hash": -871680168985845205, "size": 3, "sort": "object"}, "hash": -1146585631412276488}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/id": 10}, "hash": -6250250439178439933, "size": 3, "sort": "object"}, "hash": 2331945501239501472}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 2}, "c2": {"py/id": 10}, "hash": 6030886835341869276, "size": 3, "sort": "object"}, "hash": -4166004973050246567}]'
    return deserialize_from_string(features)


def try_clear5_features(lang):
    """ A model learnt by generating the whole 4-block state space """
    features = '[{"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -6921336085248452788, "name": "holding", "size": 1, "sort": "object"}, "hash": 3064685275814954734}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -7025971841800681580, "name": "clear", "size": 1, "sort": "object"}, "c2": {"py/object": "tarski.dl.concepts.NominalConcept", "hash": 41964751808569603, "name": "a", "size": 1, "sort": "object"}, "hash": -8101924824208142983, "size": 3, "sort": "object"}, "hash": -3172731559657890995}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/object": "tarski.dl.concepts.ForallConcept", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/id": 6}, "hash": 138301352325322456, "r": {"py/object": "tarski.dl.concepts.PrimitiveRole", "hash": -3332798531935250846, "name": "on", "size": 1, "sort": ["object", "object"]}, "size": 3, "sort": "object"}, "hash": 384461634614258817, "r": {"py/object": "tarski.dl.concepts.StarRole", "hash": 731930046869475627, "r": {"py/object": "tarski.dl.concepts.InverseRole", "hash": 2692933858444408910, "r": {"py/id": 11}, "size": 2, "sort": ["object", "object"]}, "size": 3, "sort": {"py/id": 15}}, "size": 7, "sort": "object"}, "c2": {"py/object": "tarski.dl.concepts.NotConcept", "c": {"py/id": 5}, "hash": -5348823898257502101, "size": 2, "sort": "object"}, "hash": -5596910201207631836, "size": 10, "sort": "object"}, "hash": 8757270176812848582}]'
    return deserialize_from_string(features)


if __name__ == "__main__":
    main()
