#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from sltp.incremental import IncrementalExperiment

from defaults import generate_experiment
from common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, ijcai_paper_bw_feature_namer, \
    add_bw_domain_parameters_2, build_on_x_y_feature_set, generate_features_n_ab, get_on_x_y_feature, \
    features_clear_x
from sltp.util.misc import update_dict
from sltp.util.serialization import deserialize_from_string


def experiment(experiment_name=None):
    domain_dir = "blocks"
    domain = "domain.pddl"

    exps = dict()

    # A small testing instance nonetheless producing an abstraction
    exps["debugging_test"] = dict(
        instances="instance_4_clear_x.pddl",
        num_states=20, num_sampled_states=10,
        max_concept_size=3, max_concept_grammar_iterations=1,
        parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # Learns a simple action model which is however overfit to 3 blocks,
    # and not sound in general
    exps["simple_clear_3"] = dict(
        instances="instance_3_clear_x.pddl",
        num_states=100, max_concept_size=3, max_concept_grammar_iterations=None,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["simple_clear_3_gc"] = update_dict(exps["simple_clear_3"], parameter_generator=None)
    exps["simple_clear_3_gc_blai"] = update_dict(exps["simple_clear_3_gc"], pipeline="maxsat_poly")

    # This example shows that with 4 blocks, even if we expand all states, the model is still overfit to the 4 blocks
    exps["simple_clear_4"] = dict(
        instances="instance_4_clear_x.pddl",
        num_states=150, max_concept_size=10, max_concept_grammar_iterations=3, num_sampled_states=None,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    # With these settings we generate the desired m(x):
    # card[And(And(Forall(Star(on),Not({a})), Forall(Star(Inverse(on)),Not({a}))), And(Not(holding), Not({a})))] 18
    # And indeed learnt the correct state space!!!
    exps["aaai_ijcai_features_on_clear_5_rnd"] = dict(
        instances=[
            "inst_clear_x_1.pddl",
            # "inst_clear_x_2.pddl",
        ],
        num_states=1000, num_sampled_states=100, random_seed=10,
        max_concept_size=18, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["simple_clear_3_completeness_opt"] = update_dict(exps["simple_clear_3"],
                                                  instances=["instance_3_clear_x_2.pddl"]*2,
                                                  num_states=500, max_width=[-1, 2],
                                                  num_sampled_states=100,
                                                  complete_only_wrt_optimal=True,
                                                  # enforce_features=get_on_x_y_feature
                                                  )

    exps["aaai_clear_x_simple_hybrid"] = dict(
        instances="instance_5_clear_x_1.pddl",
        test_domain=domain, test_instances=["instance_5_clear_x_2.pddl"],
        num_states=2000, max_width=[-1],
        num_sampled_states=300,
        complete_only_wrt_optimal=True,
        max_concept_size=8, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,
    )
    exps["aaai_clear_x_blai"] = update_dict(exps["aaai_clear_x_simple_hybrid"], wsat_solver_verbose=True, pipeline="maxsat_poly")
    
    # Same but using goal-concepts instead of goal parameters:
    exps["aaai_clear_x_simple_hybrid_gc"] = update_dict(exps["aaai_clear_x_simple_hybrid"], parameter_generator=None)

    exps["aaai_clear_x_8blocks"] = update_dict(exps["aaai_clear_x_simple_hybrid"],
                                       instances="instance_8_clear_x_0.pddl", )

    exps["aaai_clear_x_no_marking"] = update_dict(exps["aaai_clear_x_simple_hybrid"],
                                          complete_only_wrt_optimal=False,  # num_sampled_states=200,
                                          # concept_generator=build_ijcai_paper_bw_concepts,
                                          )

    exps["aaai_clear_x_no_marking_blai"] = update_dict(
        exps["aaai_clear_x_no_marking"], pipeline="maxsat_poly",)

    exps["aaai_clear_x_no_marking_8blocks"] = update_dict(exps["aaai_clear_x_no_marking"],
                                                  instances="instance_8_clear_x_0.pddl",)

    exps["aaai_clear_x_no_marking_k18"] = update_dict(exps["aaai_clear_x_no_marking"],
                                              complete_only_wrt_optimal=False,  # num_sampled_states=200,
                                              max_concept_size=18,
                                              )

    exps["clear_x_incremental"] = dict(
        # domain_dir="blocks-downward",
        experiment_class=IncrementalExperiment,
        # instances=["probBLOCKS-5-0.pddl", "probBLOCKS-6-0.pddl", "probBLOCKS-7-0.pddl"],
        # instances=["probBLOCKS-4-0.pddl", "probBLOCKS-4-1.pddl", "probBLOCKS-4-2.pddl"],
        instances=["instance_8_clear_x_0.pddl", "instance_5_clear_x.pddl"],  # , "instance_4_clear_x.pddl"],
        # instances=["probBLOCKS-4-0.pddl"],
        num_states=5000, max_width=-1,
        initial_sample_size=5,
        # num_sampled_states=50,
        complete_only_wrt_optimal=False,
        max_concept_grammar_iterations=3,
        initial_concept_bound=16, max_concept_bound=16, concept_bound_step=2,
        batch_refinement_size=1,
        # quiet=True,
        clean_workspace=False,
        concept_generator=None,
        parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,
    )

    # Learns a simple action model which is however overfit to 3 blocks, and not sound in general
    exps["toy_clear_incremental"] = dict(
        experiment_class=IncrementalExperiment,
        instances=["instance_3_clear_x.pddl"],
        max_concept_grammar_iterations=3,
        num_states=100,
        quiet=True,
        initial_sample_size=100,
        initial_concept_bound=10, max_concept_bound=10, concept_bound_step=2,
        concept_generator=None, parameter_generator=add_bw_domain_parameters,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["aaai_clear_x_no_marking_2"] = update_dict(
        exps["aaai_clear_x_no_marking"],
        instances=["instance_5_clear_x_1.pddl"],#,"instance_5_clear_x_2.pddl",],
        num_states=2000, max_width=[2],
        num_sampled_states=100,
    )

    #
    exps["ijcai_features_on_clear_5"] = update_dict(exps["aaai_ijcai_features_on_clear_5_rnd"], num_sampled_states=None)

    # Goal here is on(x,y).
    exps["bw_on_x_y_4"] = dict(
        instances="instance_4_on_x_y.pddl",
        num_states=200, num_sampled_states=None, random_seed=12,
        max_concept_size=31, max_concept_grammar_iterations=4,
        concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["bw_on_x_y_4_rnd"] = update_dict(exps["bw_on_x_y_4"], num_sampled_states=60)

    # Goal here is on(x,y).
    exps["bw_on_x_y_5"] = dict(
        instances="instance_5_on_x_y.pddl",
        num_states=1000, num_sampled_states=None, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["bw_on_x_y_5_iw"] = dict(
        instances=["instance_9_on_x_y_1.pddl", "instance_9_on_x_y_2.pddl", "holding_a_b_unclear.pddl"],#, "instance_9_on_x_y_4.pddl"],
        num_states=1000, max_width=2,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["aaai_bw_on_x_y_completeness_opt"] = update_dict({},
                                                  instances=[
                                                      # "inst_on_x_y_10.pddl",
                                                      # "inst_on_x_y_11.pddl",
                                                      # "inst_on_x_y_12.pddl",
                                                      # "inst_on_x_y_13.pddl",

                                                      # "inst_on_x_y_15.pddl",
                                                      "inst_on_x_y_16.pddl",
                                                      "inst_on_x_y_14.pddl",
                                                      "holding_a_b_unclear.pddl",

                                                  ],
                                                  num_states=2000, max_width=[-1],
                                                  num_sampled_states=[50, 50, 1],
                                                  complete_only_wrt_optimal=True,
                                                  sampling="all",
                                                  # enforce_features=get_on_x_y_feature,
                                                  # feature_generator=features_clear_x,
                                                  max_concept_size=8, max_concept_grammar_iterations=3,
                                                  concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                                  feature_namer=ijcai_paper_bw_feature_namer,
                                                  )

    exps["aaai_bw_on_x_y_completeness_opt_blai"] = update_dict(exps["aaai_bw_on_x_y_completeness_opt"],
                                                               pipeline="maxsat_poly")

    exps["aaai_bw_on_x_y_completeness_opt_no_marking"] = update_dict(exps["aaai_bw_on_x_y_completeness_opt"],
                                                             complete_only_wrt_optimal=False)
    
    exps["aaai_bw_on_x_y_completeness_opt_no_marking_blai"] = update_dict(
        exps["aaai_bw_on_x_y_completeness_opt_no_marking"], pipeline="maxsat_poly",)

    exps["bw_on_x_y_dt_iw"] = dict(
        instances=["on_x_y_dt_1.pddl", "holding_a_b_unclear.pddl"],
        num_states=49999, max_width=2,
        max_concept_size=10, max_concept_grammar_iterations=3,
        parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["validate_bw_on_x_y_dt_iw"] = dict(
        instances=["on_x_y_dt_1.pddl"],
        num_states=49999, max_width=2,
        max_concept_size=10, max_concept_grammar_iterations=3,
        parameter_generator=add_bw_domain_parameters_2,
        feature_generator=generate_features_n_ab,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["bw_on_x_y_dt_iw_fixed_goal"] = update_dict(exps["bw_on_x_y_dt_iw"],
                                             instances=["on_x_y_dt_1.pddl"],
                                             enforce_features=get_on_x_y_feature)

    exps["bw_on_x_y_dt_completeness_opt"] = update_dict(exps["bw_on_x_y_dt_iw"],
                                                instances=["on_x_y_dt_1.pddl"],
                                                num_states=50000, max_width=2,
                                                num_sampled_states=14000,
                                                complete_only_wrt_optimal=True,
                                                # feature_generator=weird_feature_set,
                                                enforce_features=get_on_x_y_feature)

    exps["bw_on_x_y_dt_completeness_opt2"] = update_dict(exps["bw_on_x_y_dt_completeness_opt"],
                                                 instances=["on_x_y_dt_2.pddl", "on_x_y_dt_2.pddl"],
                                                 num_states=50000, max_width=[2, -1],
                                                 # feature_generator=weird_feature_set,
                                                 # enforce_features=None,
                                                 num_sampled_states=3000)

    exps["check_ijcai_features_on_x_y"] = dict(
        instances="instance_5_on_x_y.pddl",
        num_states=1000, max_concept_size=1, max_concept_grammar_iterations=1, num_sampled_states=100, random_seed=2,
        concept_generator=build_on_x_y_feature_set, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    exps["generate_ijcai_features_on_x_y"] = dict(
        instances="instance_5_on_x_y.pddl",
        num_states=1000, max_concept_size=34, max_concept_grammar_iterations=4, num_sampled_states=50, random_seed=2,
        concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
        feature_namer=ijcai_paper_bw_feature_namer,)

    if experiment_name not in exps:
        raise RuntimeError('No experiment named "{}" in current experiment script'.format(experiment_name))
    parameters = exps[experiment_name]
    parameters["domain_dir"] = parameters.get("domain_dir", domain_dir)
    parameters["domain"] = parameters.get("domain", domain)
    return generate_experiment(**parameters)


def deserialize_features(feature_name):
    """ A model learnt by generating the whole 4-block state space """
    serialized_features = dict(
        clear4_features='[{"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -1277705213948677935, "name": "holding", "size": 1, "sort": "object"}, "hash": -2841485854126214806}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.UniversalConcept", "hash": -9223372036852258203, "size": 0, "sort": "object"}, "hash": -7885522037927662577, "r": {"py/object": "tarski.dl.concepts.PrimitiveRole", "hash": -3816393078228064225, "name": "on", "size": 1, "sort": ["object", "object"]}, "size": 2, "sort": "object"}, "hash": 911132480848590796}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NominalConcept", "hash": -7401258638559185829, "name": "a", "size": 1, "sort": "object"}, "hash": 1733562239334784851, "r": {"py/object": "tarski.dl.concepts.StarRole", "hash": 2196893821251177549, "r": {"py/id": 6}, "size": 2, "sort": {"py/id": 7}}, "size": 4, "sort": "object"}, "hash": 8416726960899686352}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/object": "tarski.dl.concepts.NotConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 5962884402008325778, "name": "ontable", "size": 1, "sort": "object"}, "hash": 2613871271492961489, "size": 2, "sort": "object"}, "hash": 3123638717406933755, "r": {"py/id": 6}, "size": 4, "sort": "object"}, "hash": -8435674864905251064}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": 3655234103021149859, "name": "clear", "size": 1, "sort": "object"}, "hash": -871680168985845205, "size": 3, "sort": "object"}, "hash": -1146585631412276488}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 15}, "c2": {"py/id": 10}, "hash": -6250250439178439933, "size": 3, "sort": "object"}, "hash": 2331945501239501472}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/id": 2}, "c2": {"py/id": 10}, "hash": 6030886835341869276, "size": 3, "sort": "object"}, "hash": -4166004973050246567}]',
        clear5_features='[{"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -6921336085248452788, "name": "holding", "size": 1, "sort": "object"}, "hash": 3064685275814954734}, {"py/object": "tarski.dl.features.EmpiricalBinaryConcept", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/object": "tarski.dl.concepts.PrimitiveConcept", "hash": -7025971841800681580, "name": "clear", "size": 1, "sort": "object"}, "c2": {"py/object": "tarski.dl.concepts.NominalConcept", "hash": 41964751808569603, "name": "a", "size": 1, "sort": "object"}, "hash": -8101924824208142983, "size": 3, "sort": "object"}, "hash": -3172731559657890995}, {"py/object": "tarski.dl.features.ConceptCardinalityFeature", "c": {"py/object": "tarski.dl.concepts.AndConcept", "c1": {"py/object": "tarski.dl.concepts.ForallConcept", "c": {"py/object": "tarski.dl.concepts.ExistsConcept", "c": {"py/id": 6}, "hash": 138301352325322456, "r": {"py/object": "tarski.dl.concepts.PrimitiveRole", "hash": -3332798531935250846, "name": "on", "size": 1, "sort": ["object", "object"]}, "size": 3, "sort": "object"}, "hash": 384461634614258817, "r": {"py/object": "tarski.dl.concepts.StarRole", "hash": 731930046869475627, "r": {"py/object": "tarski.dl.concepts.InverseRole", "hash": 2692933858444408910, "r": {"py/id": 11}, "size": 2, "sort": ["object", "object"]}, "size": 3, "sort": {"py/id": 15}}, "size": 7, "sort": "object"}, "c2": {"py/object": "tarski.dl.concepts.NotConcept", "c": {"py/id": 5}, "hash": -5348823898257502101, "size": 2, "sort": "object"}, "hash": -5596910201207631836, "size": 10, "sort": "object"}, "hash": 8757270176812848582}]',
    )
    feature_set = deserialize_from_string(serialized_features[feature_name])

    def generator(lang):
        return feature_set
    return generator


def weird_feature_set(lang):
    tower_features1 = '[{"c": {"sort": "object", "py/object": "tarski.dl.concepts.AndConcept", "c1": {"sort": "object", "size": 3, "py/object": "tarski.dl.concepts.ExistsConcept", "c": {"size": 1, "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object", "name": "b", "hash": -8502480781994927151}, "r": {"size": 1, "py/object": "tarski.dl.concepts.PrimitiveRole", "sort": ["object", "object"], "name": "on", "hash": 5776189647135905784}, "hash": 2674091572547234808}, "c2": {"size": 1, "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object", "name": "a", "hash": -2014675405248281288}, "size": 5, "hash": -9035486624804985343}, "py/object": "tarski.dl.features.EmpiricalBinaryConcept", "hash": -1364766800539080060}, {"c": {"size": 1, "py/object": "tarski.dl.concepts.PrimitiveConcept", "sort": "object", "name": "holding", "hash": 7229331559112486443}, "py/object": "tarski.dl.features.EmpiricalBinaryConcept", "hash": -540400808086029706}, {"c": {"sort": "object", "size": 8, "py/object": "tarski.dl.concepts.ExistsConcept", "c": {"sort": "object", "py/object": "tarski.dl.concepts.AndConcept", "c1": {"size": 1, "py/object": "tarski.dl.concepts.PrimitiveConcept", "sort": "object", "name": "clear", "hash": -739697115110111347}, "c2": {"c": {"size": 1, "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object", "name": "b", "hash": -8502480781994927151}, "size": 2, "py/object": "tarski.dl.concepts.NotConcept", "sort": "object", "hash": 8756838853437654867}, "size": 4, "hash": -1298651480141169035}, "r": {"size": 3, "py/object": "tarski.dl.concepts.StarRole", "sort": {"py/id": 20}, "r": {"size": 2, "py/object": "tarski.dl.concepts.InverseRole", "sort": ["object", "object"], "r": {"size": 1, "py/object": "tarski.dl.concepts.PrimitiveRole", "sort": ["object", "object"], "name": "on", "hash": 5776189647135905784}, "hash": 6322400296151783085}, "hash": 1055936545429941421}, "hash": 2197512904378287901}, "py/object": "tarski.dl.features.ConceptCardinalityFeature", "hash": -3579927205244512055}]'
    tower_features2 = '[{"hash": 6949180927697451521, "c": {"c1": {"size": 3, "r": {"size": 1, "hash": -2740289688636584686, "name": "on", "py/object": "tarski.dl.concepts.PrimitiveRole", "sort": ["object", "object"]}, "c": {"size": 1, "hash": -3218902714185107259, "name": "b", "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object"}, "sort": "object", "hash": -8028781711036654814, "py/object": "tarski.dl.concepts.ExistsConcept"}, "hash": 1384931444341254242, "sort": "object", "size": 5, "py/object": "tarski.dl.concepts.AndConcept", "c2": {"size": 1, "hash": -7897084271728384726, "name": "a", "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object"}}, "py/object": "tarski.dl.features.EmpiricalBinaryConcept"}, {"hash": -867379600521222222, "py/object": "tarski.dl.features.NullaryAtomFeature", "atom": {"hash": -2168902777557022509, "depth": 0, "name": "handempty", "py/object": "tarski.dl.concepts.NullaryAtom"}}, {"hash": 7888742611717759350, "c": {"size": 8, "r": {"py/id": 17}, "c": {"size": 6, "r": {"hash": -1329032853697446873, "r": {"hash": 1042558110122865899, "r": {"size": 1, "hash": -2740289688636584686, "name": "on", "py/object": "tarski.dl.concepts.PrimitiveRole", "sort": ["object", "object"]}, "size": 2, "py/object": "tarski.dl.concepts.InverseRole", "sort": ["object", "object"]}, "size": 3, "py/object": "tarski.dl.concepts.StarRole", "sort": {"py/id": 19}}, "c": {"size": 2, "hash": 4195320596587080666, "c": {"size": 1, "hash": -3218902714185107259, "name": "b", "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object"}, "py/object": "tarski.dl.concepts.NotConcept", "sort": "object"}, "sort": "object", "hash": -7579992454195058621, "py/object": "tarski.dl.concepts.ForallConcept"}, "sort": "object", "hash": 3102714834696527308, "py/object": "tarski.dl.concepts.ExistsConcept"}, "py/object": "tarski.dl.features.ConceptCardinalityFeature"}, {"hash": -6016550384825177845, "c": {"c1": {"size": 4, "r": {"hash": 1042414231418761324, "r": {"py/id": 17}, "size": 2, "py/object": "tarski.dl.concepts.StarRole", "sort": {"py/id": 18}}, "c": {"size": 1, "hash": -7897084271728384726, "name": "a", "py/object": "tarski.dl.concepts.NominalConcept", "sort": "object"}, "sort": "object", "hash": 5648141748725744941, "py/object": "tarski.dl.concepts.ExistsConcept"}, "hash": 8293571504602674188, "sort": "object", "size": 6, "py/object": "tarski.dl.concepts.AndConcept", "c2": {"py/id": 14}}, "py/object": "tarski.dl.features.EmpiricalBinaryConcept"}]'
    return deserialize_from_string(tower_features2)


if __name__ == "__main__":
    exp = experiment(sys.argv[1])
    exp.run(sys.argv[2:])
