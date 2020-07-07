import itertools
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
    # identified = [IdentifiedFeature(f, i, config.feature_namer(str(f)))
    #               for i, f in zip(selected_feature_ids, features)]

    policy = set()
    for (s, t) in data.sample.optimal_transitions:
        m1 = data.model_cache.get_feature_model(s)
        m2 = data.model_cache.get_feature_model(t)

        short = []
        for (i, f) in enumerate(features, start=0):
            state_val = "=0" if m1.denotation(f) == 0 else ">0"
            tx_val = str(f.diff(m1.denotation(f), m2.denotation(f))).upper()
            short += [(i, state_val), (i, tx_val)]

        policy.add(frozenset(short))

    # Minimize the DNF
    policy = minimize_dnf_policy(policy)

    # Print the policy
    logging.info("GOOD transitions:")
    for clausenum, clause in enumerate(policy, start=0):
        clause_str = []
        for fnum, val in clause:
            name = config.feature_namer(str(features[fnum]))
            clause_str.append(f'{name} {val}')

        print(f"\t{clausenum}. " + ' AND '.join(x for x in clause_str))

    return ExitCode.Success, dict(policy=policy)


def minimize_dnf_policy(dnf):
    while True:
        p1, p2, new = attempt_dnf_merge(dnf)
        if p1 is None:
            break

        # Else do the actual merge
        dnf.remove(p1)
        dnf.remove(p2)
        dnf.add(new)
    return dnf


def attempt_dnf_merge(dnf):
    for p1, p2 in itertools.combinations(dnf, 2):
        diff = p1.symmetric_difference(p2)
        diffl = list(diff)

        if len(diffl) != 2:  # More than one feature with diff value
            continue

        (f1, f1_val), (f2, f2_val) = diffl
        if f1 != f2:  # Not affecting the same feature
            continue

        if {f1_val, f2_val} == {">0", "=0"}:
            # The two conjunctions differ in that one has one literal L and the other its negation, the rest being equal
            p_merged = p1.difference(diff)
            return p1, p2, p_merged  # Meaning p1 and p2 should be merged into p_merged

    return None, None, None
