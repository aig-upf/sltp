import itertools
import logging

from .util.tools import IdentifiedFeature
from .util.tools import load_selected_features
from .returncodes import ExitCode


def compute_good_transitions(assignment, wsat_varmap_filename):
    """ Reads off the map of good-variables produced by the CNF generator, then takes the WSAT problem solution
    assignment, and computes which are the good transitions. """
    good = []
    with open(wsat_varmap_filename, 'r') as f:
        for line in f:
            var, s, t = map(int, line.rstrip().split())
            if assignment[var] is True:
                good.append((s, t))
    return list(sorted(good))


def compute_transition_separation_function(config, data, rng):
    solution = data.cnf_solution
    assert solution.solved

    # CNF variables "selected(f)" take range from 1 to num_features+1
    selected_feature_ids = [i - 1 for i in range(1, data.num_features + 1) if solution.assignment[i] is True]

    features = load_selected_features(selected_feature_ids, config.domain, config.serialized_feature_filename)
    features = [IdentifiedFeature(f, i, config.feature_namer(str(f))) for i, f in zip(selected_feature_ids, features)]

    good_transitions = compute_good_transitions(solution.assignment, config.wsat_varmap_filename)

    # debugging = []
    policy = set()
    for (s, t) in good_transitions:
        m1 = data.model_cache.get_feature_model(s)
        m2 = data.model_cache.get_feature_model(t)

        clause = []
        for f in features:
            clause += [DNFAtom(f, f.denotation(m1) != 0),
                       DNFAtom(f, f.feature.diff(f.denotation(m1), f.denotation(m2)))]

        # if ' AND '.join(sorted(map(str, clause))) == "nballs-A NILs AND nballs-A>0 AND ncarried NILs AND ncarried>0 AND nfree-grippers NILs AND nfree-grippers>0 AND robot-at-B ADDs AND robot-at-B=0":
        #     debugging.append(f"({s}, {t})")

        policy.add(frozenset(clause))

    # Minimize the DNF
    policy = minimize_dnf_policy(policy)

    # Print the policy
    # print(" ".join(debugging))
    logging.info("GOOD transitions:")
    for i, clause in enumerate(policy, start=0):
        print(f"\t{i}. " + ' AND '.join(sorted(map(str, clause))))

    return ExitCode.Success, dict(policy_dnf=TransitionSeparationPolicy(features, policy))


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

        atom1, atom2 = diffl
        if atom1.feature != atom2.feature:  # Not affecting the same feature
            continue

        if {atom1.value, atom2.value} == {True, False}:
            # The two conjunctions differ in that one has one literal L and the other its negation, the rest being equal
            p_merged = p1.difference(diff)
            return p1, p2, p_merged  # Meaning p1 and p2 should be merged into p_merged

    return None, None, None


class DNFAtom:
    def __init__(self, feature, value):
        self.feature = feature
        self.value = value

    def is_state_feature(self):
        return isinstance(self.value, bool)

    def __str__(self):
        if self.is_state_feature():
            return f'{self.feature}>0' if self.value else f'{self.feature}=0'
        # else, we have a transition feature
        return f'{self.feature} {str(self.value).upper()}s'
    __repr__ = __str__

    def __hash__(self):
        return hash((self.feature, self.value))

    def __eq__(self, other):
        return self.feature == other.feature and self.value == other.value


class TransitionSeparationPolicy:
    def __init__(self, features, policy_dnf):
        self.features = features
        self.dnf = policy_dnf

    def transition_is_good(self, m0, m1):
        for clause in self.dnf:
            all_atoms_true = True

            # If the given transition satisfies any of the clauses in the DNF, we consider it "good"
            for atom in clause:
                feat = atom.feature

                if atom.is_state_feature():
                    state_val = feat.denotation(m0) != 0
                    if state_val != atom.value:
                        all_atoms_true = False
                else:
                    tx_val = feat.feature.diff(feat.denotation(m0), feat.denotation(m1))
                    if tx_val != atom.value:
                        all_atoms_true = False

            if all_atoms_true:
                return True

        return False
