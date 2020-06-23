import logging

from .util.tools import IdentifiedFeature
from .util.tools import load_selected_features
from .returncodes import ExitCode


def compute_transition_separation_function(config, data, rng):
    solution = data.cnf_solution
    assert solution.solved

    # CNF variables "selected(f)" take range from 1 to num_features+1
    selected_feature_ids = [i - 1 for i in range(1, data.num_features + 1) if solution.assignment[i] is True]

    features = load_selected_features(selected_feature_ids, config.domain, config.serialized_feature_filename)
    identified = [IdentifiedFeature(f, i, config.feature_namer(str(f))) for i, f in zip(selected_feature_ids, features)]

    policy_dnf = set()
    for (s, t) in data.sample.optimal_transitions:
        m1 = data.model_cache.get_feature_model(s)
        m2 = data.model_cache.get_feature_model(t)
        policy_dnf.add(tuple(f"{config.feature_namer(str(f))} = {f.diff(m1.denotation(f), m2.denotation(f))}" for f in features))

    logging.info("Transitions labeled as GOOD in the computed policy:")
    for clause in policy_dnf:
        print("\t" + ' and '.join(x for x in clause))

    return ExitCode.Success, dict(policy=policy_dnf)

