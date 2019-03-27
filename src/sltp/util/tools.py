import itertools
from collections import defaultdict

from tarski.dl import ConceptCardinalityFeature, EmpiricalBinaryConcept, FeatureValueChange
from tarski.dl.features import Feature, NullaryAtomFeature

from ..language import parse_pddl
from .serialization import unserialize_features


def load_selected_features(feature_indexes, domain, filename):
    """ Load from the given filename the Feature objects with the specified indexes only """
    _, language, _ = parse_pddl(domain)
    feature_indexes = feature_indexes if isinstance(feature_indexes, set) else set(feature_indexes)
    objs = unserialize_features(language, filename, feature_indexes)
    assert len(feature_indexes) == len(objs)
    return objs


def process_features(features):
    """ Transform 'empirically binary features' in the given list into standard cardinality features.
    This is useful if we want to apply these features to unseen models, since it these models it won't
    necessarily be the case that the features are still binary. """
    return [ConceptCardinalityFeature(f.c) if isinstance(f, EmpiricalBinaryConcept) else f for f in features]


class IdentifiedFeature:
    """ A feature coupled with its unique ID and a possibly handcrafted name """
    def __init__(self, feature, id_, name):
        assert isinstance(feature, Feature)
        self.feature = feature
        self.id = id_
        self.name_ = name

    def denotation(self, model):
        return model.denotation(self.feature)
        # return self.feature.denotation(model)

    def complexity(self):
        return self.feature.complexity()

    def name(self):
        return str(self.feature) if self.name_ is None else self.name_

    def __str__(self):
        return self.name()
    __repr__ = __str__


class AbstractAction:
    def __init__(self, preconditions, effects, name=None):
        self.name = name
        self.preconditions = preconditions
        self.effects = frozenset(effects)
        self.hash = hash((self.__class__, self.preconditions, self.effects))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.preconditions == other.preconditions and self.effects == other.effects)

    def __str__(self):
        return "Action<{}; {}>".format(*prettyprint_precs_and_effects(self))

    __repr__ = __str__


def prettyprint_precs_and_effects(action):
    precs = " and ".join(sorted(print_atom(*atom) for atom in action.preconditions))
    effs = ", ".join(sorted(str(eff) for eff in action.effects))
    return precs, effs


def print_atom(feature, polarity):
    assert isinstance(feature, IdentifiedFeature)
    assert polarity in (True, False)

    if feature_type(feature) is bool:
        return feature if polarity else "not {}".format(feature)
    return "{} > 0".format(feature) if polarity else "{} = 0".format(feature)


def feature_type(feature):
    return bool if isinstance(feature, (EmpiricalBinaryConcept, NullaryAtomFeature)) else int


class Abstraction:
    def __init__(self, features, actions):
        assert all(isinstance(f, IdentifiedFeature) for f in features)
        assert all(isinstance(a, AbstractAction) for a in actions)

        self.features = features
        self.actions = actions


class ActionEffect:
    def __init__(self, feature, change):
        assert isinstance(feature, IdentifiedFeature)
        self.feature = feature
        self.change = change

    def __hash__(self):
        return hash((self.__class__, self.feature.id, self.change))

    def __eq__(self, other):
        return (self.__class__ is other.__class__ and self.feature.id == other.feature.id and
                self.change == other.change)

    def __str__(self):
        return self.print_named()

    __repr__ = __str__

    def print_named(self, namer=lambda s: s):
        name = namer(self.feature.name())
        if self.change == FeatureValueChange.ADD:
            return name
        if self.change == FeatureValueChange.DEL:
            return "NOT {}".format(name)
        if self.change == FeatureValueChange.ADD_OR_NIL:
            return "ADD* {}".format(name)
        if self.change == FeatureValueChange.INC:
            return "INC {}".format(name)
        if self.change == FeatureValueChange.DEC:
            return "DEC {}".format(name)
        if self.change == FeatureValueChange.INC_OR_NIL:
            return "INC* {}".format(name)
        raise RuntimeError("Unexpected effect type")

    def print_qnp_named(self, namer=lambda s: s):
        name = namer(self.feature.name())
        if self.change in (FeatureValueChange.ADD, FeatureValueChange.INC):
            return "{} 1".format(name)

        if self.change in (FeatureValueChange.DEL, FeatureValueChange.DEC):
            return "{} 0".format(name)

        if self.change in (FeatureValueChange.INC_OR_NIL, FeatureValueChange.ADD_OR_NIL):
            assert False, "Relaxed INC semantics not supported for QNP"
        raise RuntimeError("Unexpected effect type")


def optimize_abstract_action_model(states, actions):
    actions_grouped_by_effects = defaultdict(set)
    for a in actions:
        actions_grouped_by_effects[a.effects].add(a.preconditions)
    merged = []
    for effs, action_precs in actions_grouped_by_effects.items():
        for prec in merge_precondition_sets(action_precs):
            merged.append(AbstractAction(prec, effs))

    return states, merged


def attempt_single_merge(action_precs):
    for p1, p2 in itertools.combinations(action_precs, 2):
        diff = p1.symmetric_difference(p2)
        diffl = list(diff)
        if len(diffl) == 2 and diffl[0][0] == diffl[1][0]:
            # The two conjunctions differ in that one has one literal L and the other its negation, the rest being equal
            assert diffl[0][1] != diffl[1][1]
            p_merged = p1.difference(diff)
            return p1, p2, p_merged  # Meaning p1 and p2 should be merged into p_merged
    return None


def merge_precondition_sets(action_precs):
    action_precs = action_precs.copy()
    while True:
        res = attempt_single_merge(action_precs)
        if res is None:
            break

        # else do the actual merge
        p1, p2, new = res
        action_precs.remove(p1)
        action_precs.remove(p2)
        action_precs.add(new)
    return action_precs


def minimize_dnf(dnf):
    return merge_precondition_sets(dnf)


def generate_effect(feature, qchange):
    assert qchange in (-1, 1)

    if feature_type(feature) is bool:
        change = {-1: FeatureValueChange.DEL, 1: FeatureValueChange.ADD}[qchange]
    else:  # else type must be int
        change = {-1: FeatureValueChange.DEC, 1: FeatureValueChange.INC}[qchange]

    return ActionEffect(feature, change)
