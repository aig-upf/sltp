import itertools

from tarski.dl import FeatureValueChange

from .language import parse_pddl
from .util.serialization import unserialize_feature
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


def compute_transition_classification_policy(config, data, rng):
    if config.transition_classification_policy is not None:
        rules = config.transition_classification_policy()
        _, language, _ = parse_pddl(config.domain)
        policy = TransitionClassificationPolicy.parse(rules, language, config.feature_namer)
        return ExitCode.Success, dict(transition_classification_policy=policy)

    solution = data.cnf_solution
    assert solution.solved

    # CNF variables "selected(f)" take range from 1 to num_features+1
    selected_feature_ids = [i - 1 for i in range(1, data.num_features + 1) if solution.assignment[i] is True]

    features = load_selected_features(selected_feature_ids, config.domain, config.serialized_feature_filename)
    features = [IdentifiedFeature(f, i, config.feature_namer(str(f))) for i, f in zip(selected_feature_ids, features)]

    good_transitions = compute_good_transitions(solution.assignment, config.wsat_varmap_filename)

    policy = TransitionClassificationPolicy(features)

    for (s, t) in good_transitions:
        m1 = data.model_cache.get_feature_model(s)
        m2 = data.model_cache.get_feature_model(t)

        clause = []
        for f in features:
            clause += [DNFAtom(f, f.denotation(m1) != 0),
                       DNFAtom(f, f.feature.diff(f.denotation(m1), f.denotation(m2)))]

        policy.add_clause(frozenset(clause))

    policy.minimize()
    policy.print()
    return ExitCode.Success, dict(transition_classification_policy=policy)


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


class TransitionClassificationPolicy:
    def __init__(self, features):
        self.features = features
        self.dnf = set()

    def add_clause(self, clause):
        self.dnf.add(clause)

    def minimize(self):
        self.dnf = minimize_dnf_policy(self.dnf)

    def transition_is_good(self, m0, m1):
        # If the given transition satisfies any of the clauses in the DNF, we consider it "good"
        return any(self.does_transition_satisfy_clause(clause, m0, m1) for clause in self.dnf)

    @staticmethod
    def does_transition_satisfy_clause(clause, m0, m1):
        for atom in clause:
            feat = atom.feature

            if atom.is_state_feature():
                state_val = feat.denotation(m0) != 0
                if state_val != atom.value:
                    return False
            else:
                tx_val = feat.feature.diff(feat.denotation(m0), feat.denotation(m1))
                if tx_val != atom.value:
                    return False
        return True

    def print(self):
        print("Transition-classification policy with the following transitions labeled as good:")
        for i, clause in enumerate(self.dnf, start=0):
            print(f"  {i}. " + ' AND '.join(sorted(map(str, clause))))

    @staticmethod
    def parse(rules, language, feature_namer):
        """ Create a classification policy from a set of strings representing the clauses """
        policy = TransitionClassificationPolicy(features=[])

        allfeatures = dict()

        for clause in rules:
            atoms = []
            for feature_str, value in clause:
                f = allfeatures.get(feature_str)
                if f is None:
                    f = unserialize_feature(language, feature_str)
                    allfeatures[feature_str] = f = IdentifiedFeature(f, len(allfeatures), feature_namer(str(f)))

                # Convert the value to an object
                value = {
                    "=0": False,
                    ">0": True,
                    "INC": FeatureValueChange.INC,
                    "NIL": FeatureValueChange.NIL,
                    "DEC": FeatureValueChange.DEC,
                    "ADD": FeatureValueChange.ADD,
                    "DEL": FeatureValueChange.DEL,
                }[value]

                atoms.append(DNFAtom(f, value))

            policy.add_clause(frozenset(atoms))

        return policy
