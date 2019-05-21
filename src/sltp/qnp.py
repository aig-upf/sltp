

def generate_encoding(config, data, rng):
    config.qnp_abstraction_filename
    config.qnp_prefix

    states, actions, features = data.cnf_translator.compute_action_model(data.cnf_solution.assignment, config)
    data.cnf_translator.compute_qnp(actions, features, config, data.sample)
    # return dict(abstract_states=states, abstract_actions=actions)
    return dict()
