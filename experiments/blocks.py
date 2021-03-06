
from sltp.incremental import IncrementalExperiment

from common import build_ijcai_paper_bw_concepts, add_bw_domain_parameters, \
    add_bw_domain_parameters_2, build_on_x_y_feature_set, generate_features_n_ab, get_on_x_y_feature, \
    features_clear_x
from sltp.util.misc import update_dict
from sltp.util.names import blocksworld_names


def experiments():
    domain = "domain.pddl"
    base = dict(
        domain_dir="blocks",
        domain=domain,
    )

    exps = dict()

    # A small testing instance nonetheless producing an abstraction
    exps["debugging_test"] = update_dict(base,
                                         instances="instance_4_clear_x.pddl",
                                         num_states=20, num_sampled_states=10,
                                         max_concept_size=3, max_concept_grammar_iterations=1,
                                         parameter_generator=add_bw_domain_parameters,
                                         feature_namer=blocksworld_names, )

    # Learns a simple action model which is however overfit to 3 blocks,
    # and not sound in general
    exps["simple_clear_3"] = update_dict(base,
                                         instances="instance_3_clear_x.pddl",
                                         test_domain=domain, test_instances=["instance_8_clear_x_0.pddl"],
                                         num_states=100, max_concept_size=6, max_concept_grammar_iterations=None,
                                         concept_generator=None, parameter_generator=add_bw_domain_parameters,
                                         feature_namer=blocksworld_names, )

    exps["simple_clear_3_gc"] = update_dict(exps["simple_clear_3"], parameter_generator=None)
    exps["simple_clear_3_gc_blai"] = update_dict(exps["simple_clear_3_gc"], pipeline="maxsat_poly")

    # This example shows that with 4 blocks, even if we expand all states, the model is still overfit to the 4 blocks
    exps["simple_clear_4"] = update_dict(base,
                                         instances="instance_4_clear_x.pddl",
                                         num_states=150, max_concept_size=10, max_concept_grammar_iterations=3, num_sampled_states=None,
                                         concept_generator=None, parameter_generator=add_bw_domain_parameters,
                                         feature_namer=blocksworld_names, )

    # With these settings we generate the desired m(x):
    # card[And(And(Forall(Star(on),Not({a})), Forall(Star(Inverse(on)),Not({a}))), And(Not(holding), Not({a})))] 18
    # And indeed learnt the correct state space!!!
    exps["aaai_ijcai_features_on_clear_5_rnd"] = update_dict(base,
                                                             instances=[
            "inst_clear_x_1.pddl",
            # "inst_clear_x_2.pddl",
        ],
                                                             num_states=1000, num_sampled_states=100, random_seed=10,
                                                             max_concept_size=18, max_concept_grammar_iterations=3,
                                                             concept_generator=None, parameter_generator=add_bw_domain_parameters,
                                                             feature_namer=blocksworld_names, )

    exps["simple_clear_3_completeness_opt"] = update_dict(exps["simple_clear_3"],
                                                  instances=["instance_3_clear_x_2.pddl"]*2,
                                                  num_states=500, max_width=[-1, 2],
                                                  num_sampled_states=100,
                                                  complete_only_wrt_optimal=True,
                                                  # enforce_features=get_on_x_y_feature
                                                  )

    exps["aaai_clear_x_simple_hybrid"] = update_dict(base,
                                                     instances="instance_5_clear_x_1.pddl",
                                                     test_domain=domain, test_instances=["instance_5_clear_x_2.pddl"],
                                                     num_states=2000, max_width=[-1],
                                                     num_sampled_states=300,
                                                     complete_only_wrt_optimal=True,
                                                     max_concept_size=8, max_concept_grammar_iterations=3,
                                                     concept_generator=None, parameter_generator=add_bw_domain_parameters,
                                                     feature_namer=blocksworld_names,
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

    exps["clear_x_incremental"] = update_dict(base,
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
                                              feature_namer=blocksworld_names,
                                              )

    # Learns a simple action model which is however overfit to 3 blocks, and not sound in general
    exps["toy_clear_incremental"] = update_dict(base,
                                                experiment_class=IncrementalExperiment,
                                                instances=["instance_3_clear_x.pddl"],
                                                max_concept_grammar_iterations=3,
                                                num_states=100,
                                                quiet=True,
                                                initial_sample_size=100,
                                                initial_concept_bound=10, max_concept_bound=10, concept_bound_step=2,
                                                concept_generator=None, parameter_generator=add_bw_domain_parameters,
                                                feature_namer=blocksworld_names, )

    exps["aaai_clear_x_no_marking_2"] = update_dict(
        exps["aaai_clear_x_no_marking"],
        instances=["instance_5_clear_x_1.pddl"],#,"instance_5_clear_x_2.pddl",],
        num_states=2000, max_width=[2],
        num_sampled_states=100,
    )

    #
    exps["ijcai_features_on_clear_5"] = update_dict(exps["aaai_ijcai_features_on_clear_5_rnd"], num_sampled_states=None)

    # Goal here is on(x,y).
    exps["bw_on_x_y_4"] = update_dict(base,
                                      instances="instance_4_on_x_y.pddl",
                                      num_states=200, num_sampled_states=None, random_seed=12,
                                      max_concept_size=31, max_concept_grammar_iterations=4,
                                      concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                      feature_namer=blocksworld_names, )

    exps["bw_on_x_y_4_rnd"] = update_dict(exps["bw_on_x_y_4"], num_sampled_states=60)

    # Goal here is on(x,y).
    exps["bw_on_x_y_5"] = update_dict(base,
                                      instances="instance_5_on_x_y.pddl",
                                      num_states=1000, num_sampled_states=None, random_seed=12,
                                      max_concept_size=10, max_concept_grammar_iterations=3,
                                      concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                      feature_namer=blocksworld_names, )

    exps["bw_on_x_y_5_gc"] = update_dict(exps["bw_on_x_y_5"], parameter_generator=None,)

    exps["bw_on_x_y_5_iw"] = update_dict(base,
                                         instances=["instance_9_on_x_y_1.pddl", "instance_9_on_x_y_2.pddl", "holding_a_b_unclear.pddl"],  #, "instance_9_on_x_y_4.pddl"],
                                         num_states=1000, max_width=2,
                                         max_concept_size=10, max_concept_grammar_iterations=3,
                                         concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                         feature_namer=blocksworld_names, )

    exps["aaai_bw_on_x_y_completeness_opt"] = update_dict(base,
                                                          instances=[
                                                      "inst_on_x_y_16.pddl",
                                                      "inst_on_x_y_14.pddl",
                                                      "holding_a_b_unclear.pddl",

                                                  ],
                                                          # test_domain=domain,
                                                          # We cannot test this, since works only for states where A and B are on diff towers
                                                          # test_instances=["inst_on_x_y_7.pddl"],
                                                          num_states=2000, max_width=[-1],
                                                          num_sampled_states=[50, 50, 1],
                                                          complete_only_wrt_optimal=True,
                                                          sampling="all",
                                                          # enforce_features=get_on_x_y_feature,
                                                          # feature_generator=features_clear_x,
                                                          max_concept_size=8, max_concept_grammar_iterations=3,
                                                          concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                                          feature_namer=blocksworld_names,
                                                          )

    exps["aaai_bw_on_x_y_completeness_opt_blai"] = update_dict(exps["aaai_bw_on_x_y_completeness_opt"],
                                                               pipeline="maxsat_poly")

    exps["aaai_bw_on_x_y_completeness_opt_no_marking"] = update_dict(exps["aaai_bw_on_x_y_completeness_opt"],
                                                             complete_only_wrt_optimal=False)
    
    exps["aaai_bw_on_x_y_completeness_opt_no_marking_blai"] = update_dict(
        exps["aaai_bw_on_x_y_completeness_opt_no_marking"], pipeline="maxsat_poly",)

    exps["bw_on_x_y_dt_iw"] = update_dict(base,
                                          instances=["on_x_y_dt_1.pddl", "holding_a_b_unclear.pddl"],
                                          num_states=49999, max_width=2,
                                          max_concept_size=10, max_concept_grammar_iterations=3,
                                          parameter_generator=add_bw_domain_parameters_2,
                                          feature_namer=blocksworld_names, )

    exps["validate_bw_on_x_y_dt_iw"] = update_dict(base,
                                                   instances=["on_x_y_dt_1.pddl"],
                                                   num_states=49999, max_width=2,
                                                   max_concept_size=10, max_concept_grammar_iterations=3,
                                                   parameter_generator=add_bw_domain_parameters_2,
                                                   feature_generator=generate_features_n_ab,
                                                   feature_namer=blocksworld_names, )

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

    exps["check_ijcai_features_on_x_y"] = update_dict(base,
                                                      instances="instance_5_on_x_y.pddl",
                                                      num_states=1000, max_concept_size=1, max_concept_grammar_iterations=1, num_sampled_states=100, random_seed=2,
                                                      concept_generator=build_on_x_y_feature_set, parameter_generator=add_bw_domain_parameters_2,
                                                      feature_namer=blocksworld_names, )

    exps["generate_ijcai_features_on_x_y"] = update_dict(base,
                                                         instances="instance_5_on_x_y.pddl",
                                                         num_states=1000, max_concept_size=34, max_concept_grammar_iterations=4, num_sampled_states=50, random_seed=2,
                                                         concept_generator=None, parameter_generator=add_bw_domain_parameters_2,
                                                         feature_namer=blocksworld_names, )

    exps["clear_two_atoms"] = update_dict(base,
                                          instances="instance_5_clear_x_y.pddl",
                                          test_domain=domain, test_instances=["instance_5_clear_x_2.pddl"],
                                          num_states=2000, max_width=[-1],
                                          num_sampled_states=300,
                                          complete_only_wrt_optimal=True,
                                          max_concept_size=8,
                                          concept_generator=None,
                                          parameter_generator=None,
                                          feature_namer=blocksworld_names,
                                          )

    return exps
