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

    # Learns a simple action model which is however overfit to 3 blocks,
    # and not sound in general
    simple_clear_3 = dict(
        instance="instance_3_clear_x.pddl",
        num_states=100, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=concept_generator, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    #
    simple_clear_4 = dict(
        instance="instance_4_clear_x.pddl",
        num_states=150, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_generator=try_clear4_features,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # With these settings we generate the desired m(x):
    # card[And(And(Forall(Star(on),Not({a})), Forall(Star(Inverse(on)),Not({a}))), And(Not(holding), Not({a})))] 18
    can_we_learn_ijcai_features_on_clear_5 = dict(
        instance="instance_5_clear_x_1.pddl",
        num_states=70, max_concept_size=18, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    check_ijcai_features_on_clear_5 = dict(
        instance="instance_5_clear_x.pddl",
        num_states=600, max_concept_size=1, max_concept_grammar_iterations=1,
        concept_generator=build_ijcai_paper_bw_concepts, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exp = generate_experiment(domain_dir, domain, **can_we_learn_ijcai_features_on_clear_5)
    exp.run()


def try_clear4_features(lang):
    features = '[{"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -1277705213948677935, "name": "holding", "size": 1, "sort": "object"}, "hash": -2841485854126214806}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.UniversalConcept", "hash": -9223372036852258203, "size": 0, "sort": "object"}, "hash": -7885522037927662577, "r": {"py/object": "tarski.dl.concepts.PrimitiveRole", "hash": -3816393078228064225, "name": "on", "size": 1, "sort": ["object", "object"]}, "size": 2, "sort": "object"}, "hash": 911132480848590796}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NominalConcept", "hash": -7401258638559185829, "name": "a", "size": 1, "sort": "object"}, "hash": 1733562239334784851, "r": {"py/object": "tarski.dl.concepts.StarRole", "hash": 2196893821251177549, "r": {"py/id": 6}, "size": 2, "sort": {"py/id": 7}}, "size": 4, "sort": "object"}, "hash": 8416726960899686352}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NotConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 5962884402008325778, "name": "ontable", "size": 1, "sort": "object"}, "hash": 2613871271492961489, "size": 2, "sort": "object"}, "hash": 3123638717406933755, "r": {"py/id": 6}, "size": 4, "sort": "object"}, "hash": -8435674864905251064}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 3655234103021149859, "name": "clear", "size": 1, "sort": "object"}, "hash": -871680168985845205, "size": 3, "sort": "object"}, "hash": -1146585631412276488}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/id": 10}, "hash": -6250250439178439933, "size": 3, "sort": "object"}, "hash": 2331945501239501472}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 2}, "c2": {"py/id": 10}, "hash": 6030886835341869276, "size": 3, "sort": "object"}, "hash": -4166004973050246567}]'
    return deserialize_from_string(features)


if __name__ == "__main__":
    main()
