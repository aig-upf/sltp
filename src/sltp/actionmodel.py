
from .learn_actions import create_maxsat_translator
from .returncodes import ExitCode


def compute_action_model(config, data, rng):
    assert data.cnf_solution.solved
    selected_features = data.cnf_translator.decode_solution(data.cnf_solution.assignment)
    states, actions, features = data.cnf_translator.compute_action_model(selected_features, config)
    data.cnf_translator.compute_qnp(states, actions, features, config, data)
    return ExitCode.Success, dict(abstract_actions=actions, selected_features=features)


def compute_action_model_from_feature_idxs(config, data, rng):
    sat_feature_mapping = data.sat_feature_mapping
    # Remap to the old indexes, which are the ones used in the Features vector
    selected_features = sorted(sat_feature_mapping[i] for i in data.selected_feature_idxs)
    # We create the translator just to compute the action model. TODO It'd be better to fully decouple both steps,
    # so that we don't need to perform unnecessary initialization operations here.
    translator = create_maxsat_translator(config, data)
    states, actions, features = translator.compute_action_model(selected_features, config)
    # data.cnf_translator.compute_qnp(states, actions, features, config, data)
    return ExitCode.Success, dict(abstract_actions=actions, selected_features=features)