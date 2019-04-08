import logging
import os

import numpy as np

from .util.tools import IdentifiedFeature, AbstractAction, Abstraction, optimize_abstract_action_model, \
    prettyprint_precs_and_effects, generate_effect
from .errors import CriticalPipelineError
from .sampling import TransitionSample
from .learn_actions import create_maxsat_translator, compute_completeness_info
from .util.tools import load_selected_features
from .returncodes import ExitCode
from .validator import AbstractionValidator


def compute_action_model(config, data, rng):
    solution = data.cnf_solution
    assert solution.solved
    selected_feature_ids = data.cnf_translator.decode_solution(solution.assignment)

    # Compute information relevant to the type of abstraction we want
    cinfo = compute_completeness_info(data.sample, config.complete_only_wrt_optimal)

    return _compute_abstract_action_model(config, data, selected_feature_ids, cinfo)


def compute_action_model_from_feature_idxs(config, data, rng):
    sat_feature_mapping = data.sat_feature_mapping
    # Remap to the old indexes, which are the ones used in the Features vector
    selected_features = sorted(sat_feature_mapping[i] for i in data.selected_feature_idxs)
    # We create the translator just to compute the action model. TODO It'd be better to fully decouple both steps,
    # so that we don't need to perform unnecessary initialization operations here.
    translator, sample = create_maxsat_translator(config, data.sample)
    states, actions, features = translator.compute_action_model(selected_features, config)
    # data.cnf_translator.compute_qnp(states, actions, features, config, data)
    return ExitCode.Success, dict(abstract_actions=actions, selected_features=features, post_cnf_sample=sample)


def _compute_policy_from_abstraction(config, data, abstraction, cinfo):
    logging.info("Computing policy from abstraction")
    sample = data.sample
    features = abstraction.features
    policy = dict()
    validator = AbstractionValidator(data.model_cache, sample, None)
    feature_idx = validator.compute_feature_idx(abstraction.actions)
    for s, sprime in cinfo.optimal_transitions:
        model = data.model_cache.get_feature_model(s)
        cache = {}
        app = [i for i, a in enumerate(abstraction.actions, start=1) if validator.is_applicable(cache, s, a) and
               validator.action_captures(cache, s, sprime, feature_idx[a], features)]

        if not app:
            raise RuntimeError("Training set transition not captures by any action in the abstraction")

        if len(app) > 1:
            logging.warning("Two abstract actions applicable in same abstract state")  # This shouldn't happen

        else:
            abstract_s = abstract_state(model, features)
            right_action = app[0]
            previous = policy.get(abstract_s, None)
            if previous is not None and previous != right_action:
                logging.warning("Two actions appear optimal in data-driven policy: {} and {}".
                                format(previous, right_action))  # This could happen
            policy[abstract_s] = right_action

    logging.info("Computed policy of size {}".format(len(policy)))
    print("\t" + "\n\t".join("{}: {}".format(s, a) for s, a in policy.items()))
    return policy


def compute_abstract_action_model(config, data, rng):
    solution = data.cnf_solution
    assert solution.solved

    # CNF variables "selected(f)" take range from 1 to num_features+1
    selected_feature_ids = [i - 1 for i in range(1, data.num_features + 1) if solution.assignment[i] is True]

    # Compute information relevant to the type of abstraction we want
    cinfo = compute_completeness_info(data.sample, config.complete_only_wrt_optimal)

    result, resdata = _compute_abstract_action_model(config, data, selected_feature_ids, cinfo)
    if result != ExitCode.Success:
        resdata["policy"] = None
    else:
        resdata["policy"] = _compute_policy_from_abstraction(config, data, resdata["abstraction"], cinfo)
    return result, resdata


def _compute_abstract_action_model(config, data, feature_ids, cinfo):
    features = load_selected_features(feature_ids, config.domain, config.serialized_feature_filename)
    identified = [IdentifiedFeature(f, i, config.feature_namer(str(f))) for i, f in zip(feature_ids, features)]

    # Compute the abstract state space and store it to disk
    states, actions = compute_abstraction(data.sample, cinfo, identified, data.model_cache)
    print_actions(actions, os.path.join(config.experiment_dir, 'actions.txt'))

    # Minimize the abstract action model and store it to disk
    states, actions = optimize_abstract_action_model(states, actions)
    opt_filename = os.path.join(config.experiment_dir, 'optimized.txt')
    print_actions(actions, opt_filename)
    logging.info("Minimized action model with {} actions saved in {}".format(len(actions), opt_filename))

    abstraction = Abstraction(features=identified, actions=actions)
    return ExitCode.Success, dict(abstraction=abstraction)


def compute_abstraction(sample, cinfo, features, model_cache):
    if not features:
        raise CriticalPipelineError("0-cost CNF solution - the encoding has likely some error")

    print("Features (total complexity: {}): ".format(sum(f.complexity() for f in features)))
    print('\t' + '\n\t'.join("{}. {} [k={}, id={}]".format(i, f, f.complexity(), f.id) for i, f in enumerate(features)))

    abstract_states, state_abstraction = compute_state_abstraction(sample, features, model_cache)

    abstract_actions = set()
    for s, sprime in cinfo.optimal_transitions:
        abstr_s = state_abstraction[s]
        abstr_sprime = state_abstraction[sprime]

        s_m = model_cache.get_feature_model(s)
        sprime_m = model_cache.get_feature_model(sprime)
        qchanges = [np.sign(sprime_m.denotation(f) - s_m.denotation(f)) for f in features]

        # Generate the action effects from those changes which are not NIL
        abstract_effects = [generate_effect(f, c) for f, c in zip(features, qchanges) if c != 0]

        precondition_bitmap = frozenset(zip(features, abstr_s))
        abstract_actions.add(AbstractAction(precondition_bitmap, abstract_effects))
        if len(abstract_effects) == 0:
            logging.warning("Abstract noop needed! [concrete: ({}, {}), abstract: ({}, {})]".
                            format(s, sprime, abstr_s, abstr_sprime))

    logging.info("Abstract state space: {} states and {} actions".format(len(abstract_states), len(abstract_actions)))
    return abstract_states, abstract_actions


def compute_state_abstraction(sample: TransitionSample, features, model_cache):
    abstract_states = set()
    state_abstraction = dict()
    for sid, state in sample.states.items():
        model = model_cache.get_feature_model(sid)
        abstract = abstract_state(model, features)
        abstract_states.add(abstract)
        state_abstraction[sid] = abstract
    return abstract_states, state_abstraction


def abstract_state(model, features):
    # The abstraction uses the boolean values: either F=0 or F>0
    return tuple(bool(model.denotation(f)) for f in features)


def print_actions(actions, filename):
    with open(filename, 'w') as f:
        for i, action in enumerate(actions, 1):
            precs, effs = prettyprint_precs_and_effects(action)
            action_str = "\tPRE: {}\n\tEFFS: {}".format(precs, effs)
            print("\nAction {}:\n{}".format(i, action_str), file=f)

