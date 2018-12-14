"""
Concept and feature models and related classes
"""
import copy
from bitarray import bitarray

from sltp.extensions import uncompress_extension, compress_extension
from sltp.util.misc import try_number
from tarski import Predicate, Function
from tarski.dl import UniversalConcept, EmptyConcept, NullaryAtom, PrimitiveConcept, PrimitiveRole, GoalRole, \
    GoalConcept, GoalNullaryAtom, NominalConcept, EmpiricalBinaryConcept

_STANDARD_DL_MAPPING = {0: NullaryAtom, 1: PrimitiveConcept, 2: PrimitiveRole}
_GOAL_DL_MAPPING = {0: GoalNullaryAtom, 1: GoalConcept, 2: GoalRole}
_TOP = UniversalConcept('object')
_BOT = EmptyConcept('object')


# TODO Should take into account distinction between static / dynamic facts
class DLModel:
    """ """
    def __init__(self, primitive_denotations, cache=None):
        self.primitive_denotations = primitive_denotations
        self.cache = cache
        self._universe = self.primitive_denotations[_TOP]

    def universe(self):
        return self._universe

    def denotation(self, term):
        return term.denotation(self)

    def primitive_denotation(self, term):
        return self.primitive_denotations[term]

    def compressed(self, data, arity):
        if isinstance(data, bitarray):
            return data
        return compress_extension(data, len(self._universe), arity)

    def uncompressed(self, data, arity):
        if not isinstance(data, bitarray):
            return data
        return uncompress_extension(data, len(self._universe), arity)

    def uncompressed_denotation(self, term):
        return self.uncompressed(self.denotation(term), term.ARITY)

    def compressed_denotation(self, term):
        return self.compressed(self.denotation(term), term.ARITY)


def default_denotation_by_arity(arity):  # The default values for each possible arity
    if arity == 0:
        return False
    return set()


# TODO Should take into account distinction between static / dynamic facts
class DLModelFactory:
    """ A factory of DL models tailored to a concrete universe of discourse.
    This means that probable the factory is suitable to generate models that work for a single planning instance,
    unless all of your instances happen to have the same universe. Even it that case, in the future we might want to
    implement some precomputations based on the goal and static information of the instance.
    """
    def __init__(self, universe, vocabulary, nominals, goal_conjunction=None):
        """ `vocabulary` should contain a mapping from (relevant) symbol names to the actual Tarski symbol object """
        self.universe = universe
        self.vocabulary = vocabulary  # vocabulary contains pred./function symbols relevant for any interpretation
        self.base_denotations = self.compute_base_denotations(universe, vocabulary, nominals, goal_conjunction)

    def compute_base_denotations(self, universe, vocabulary, nominals, goal_conjunction=None):
        """ Initialize the data structure that we will use to store the denotation that corresponds
            to each logical symbol in the language. For instance, the denotation of a unary symbol such as "clear"
            will be represented as a set of objects; here we just create a mapping from the "clear" symbol to
            an empty set. This incidentally enforces the closed-world assumption:
            predicates not appearing on the state trace will be assumed to have empty denotation"""
        denotations = {_TOP: universe.as_extension(), _BOT: set()}
        for p in vocabulary.values():
            ar = p.uniform_arity()
            dl_element = _STANDARD_DL_MAPPING[ar](p)
            denotations[dl_element] = default_denotation_by_arity(ar)
            if goal_conjunction is not None:
                dl_element = _GOAL_DL_MAPPING[ar](p)
                denotations[dl_element] = default_denotation_by_arity(ar)

        #
        for nominal in nominals:
            denotations[NominalConcept(nominal.symbol, nominal.sort)] = {self.universe.index(nominal.symbol)}

        # If a goal was passed, add the goal denotations (e.g. clear_G) computed from the goal "partial state"
        if goal_conjunction is not None:
            _ = [self.process_atom(atom, denotations, _GOAL_DL_MAPPING) for atom in goal_conjunction]

        return denotations

    def base_model(self):
        # TODO Might want not to create n copies of the universe?
        def process(k, v):
            if isinstance(v, bool):
                return v
            return v.copy()
        return {k: process(k, v) for k, v in self.base_denotations.items()}

    def create_model(self, state):
        """ Create a model capable of interpreting any DL concept / role under the given state """
        # Start with a copy of all the precomputed data
        denotations = self.base_model()
        _ = [self.process_atom(atom, denotations, _STANDARD_DL_MAPPING) for atom in state]
        return DLModel(denotations)

    def process_atom(self, atom, denotations, dl_mapping):
        """ Process an atom represented in format e.g. ("on", "a", "b") and add the corresponding modifications
        to the `denotations` dictionary.
        """
        assert len(atom) <= 3, "Cannot deal with arity>2 predicates or arity>1 functions yet"
        predfun = self.vocabulary[atom[0]]
        assert isinstance(predfun, (Predicate, Function))
        assert predfun.uniform_arity() == len(atom) - 1
        dl_element = dl_mapping[predfun.uniform_arity()](predfun)

        if len(atom) == 1:  # i.e. a nullary predicate
            denotations[dl_element] = True

        elif len(atom) == 2:  # i.e. a unary predicate or nullary function
            denotations[dl_element].add(self.universe.index(try_number(atom[1])))

        else:  # i.e. a binary predicate or unary function
            t = (self.universe.index(try_number(atom[1])), self.universe.index(try_number(atom[2])))
            denotations[dl_element].add(t)
            

class FeatureModel:
    """ """
    def __init__(self, concept_model):
        assert isinstance(concept_model, DLModel)
        self.concept_model = concept_model

    def denotation(self, feature):
        # TODO ADD CACHE?
        val = feature.denotation(self.concept_model)
        return bool(val) if isinstance(feature, EmpiricalBinaryConcept) else val


def compute_feature_denotation(feature, model):
    cast = bool if isinstance(feature, EmpiricalBinaryConcept) else lambda x: x
    return cast(feature.denotation(model))
