from enum import Enum

import tarski as tsk

from util.algorithms import transitive_closure, compute_min_distance


# abstract classes for concepts and roles
class Concept(object):
    ARITY = 1

    def __init__(self, sort, depth):
        self.sort = sort
        self.depth = depth

    def extension(self, cache, state, parameter_subst):
        raise NotImplementedError()

    def flatten(self):
        raise NotImplementedError()


class Role(object):
    ARITY = 2

    def __init__(self, sort, depth):
        self.sort = sort
        self.depth = depth

    def extension(self, cache, state, substitution):
        raise NotImplementedError()

    def flatten(self):
        raise NotImplementedError()


class Atom:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


class UniversalConcept(Concept):
    def __init__(self, universal_sort):
        Concept.__init__(self, universal_sort, 0)
        self.hash = hash(self.__class__)

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return self.__class__ is other.__class__

    def extension(self, cache, state, substitution):
        return cache.as_bitarray(self, state)

    def __repr__(self):
        return '<universe>'

    __str__ = __repr__

    def flatten(self):
        return [self]


class EmptyConcept(Concept):
    def __init__(self, universal_sort):
        Concept.__init__(self, universal_sort, 0)
        self.hash = hash(self.__class__)

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return self.__class__ is other.__class__

    def extension(self, cache, state, substitution):
        return cache.as_bitarray(self, state)

    def __repr__(self):
        return '<empty>'

    __str__ = __repr__

    def flatten(self):
        return [self]


class SingletonConcept(Concept):
    def __init__(self, name, sort):
        Concept.__init__(self, sort, 0)
        self.name = name
        self.hash = hash((self.__class__, self.name))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.name == other.name)

    def extension(self, cache, state, substitution):
        singleton = {cache.universe.index(self.name)}
        return cache.compress(singleton, self.ARITY)

    def __repr__(self):
        return "{{{}}}".format(self.name)

    __str__ = __repr__

    def flatten(self):
        return [self]


class BasicConcept(Concept):
    def __init__(self, predicate):
        assert isinstance(predicate, tsk.syntax.Predicate)
        Concept.__init__(self, predicate.sort[0], 0)
        self.predicate = predicate
        # This is a bit aggressive, but we assume that predicate names are unique
        self.hash = hash((self.__class__, self.name))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.name == other.name)

    @property
    def name(self):
        return self.predicate.symbol

    def extension(self, cache, state, substitution):
        return cache.as_bitarray(self, state)

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__

    def flatten(self):
        return [self]


class NotConcept(Concept):
    def __init__(self, c, universal_sort):
        assert isinstance(c, Concept)
        Concept.__init__(self, universal_sort, 1 + c.depth)
        self.c = c
        self.hash = hash((self.__class__, self.c))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.c == other.c)

    def extension(self, cache, state, substitution):
        ext_c = cache.as_bitarray(self.c, state)
        return ~ext_c

    def __repr__(self):
        return 'Not(%s)' % repr(self.c)

    __str__ = __repr__

    def flatten(self):
        return [self] + self.c.flatten()


class AndConcept(Concept):
    def __init__(self, c1, c2):
        assert isinstance(c1, Concept)
        assert isinstance(c2, Concept)
        lang = c1.sort.language
        Concept.__init__(self, most_restricted_type(lang, c1.sort, c2.sort), 1 + c1.depth + c2.depth)
        self.c1 = c1
        self.c2 = c2
        self.hash = hash((self.__class__, self.c1, self.c2))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.c1 == other.c1 and
                self.c2 == other.c2)

    def extension(self, cache, state, substitution):
        ext_c1 = cache.as_bitarray(self.c1, state)
        ext_c2 = cache.as_bitarray(self.c2, state)
        return ext_c1 & ext_c2

    def __repr__(self):
        return 'And(%s,%s)' % (repr(self.c1), repr(self.c2))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.c1.flatten() + self.c2.flatten()


class ExistsConcept(Concept):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        # The sort of an exists-concept is that of the first element of the relation
        Concept.__init__(self, r.sort[0], 1 + r.depth + c.depth)
        self.r = r
        self.c = c
        self.hash = hash((self.__class__, self.r, self.c))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.c == other.c and
                self.r == other.r)

    def extension(self, cache, state, substitution):
        ext_c = cache.as_set(self.c, state)
        ext_r = cache.as_set(self.r, state)
        # result = [x for x in objects if [z for (y, z) in ext_r if y == x and z in ext_c]]
        result = set(x for x, y in ext_r if y in ext_c)
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Exists(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r.flatten() + self.c.flatten()


class ForallConcept(Concept):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        # The sort of a forall-concept is that of the first element of the relation # TODO Check this
        Concept.__init__(self, r.sort[0], 1 + r.depth + c.depth)
        self.r = r
        self.c = c
        self.hash = hash((self.__class__, self.r, self.c))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.c == other.c and
                self.r == other.r)

    def extension(self, cache, state, substitution):
        universe = cache.universe_set
        ext_c = cache.as_set(self.c, state)
        ext_r = cache.as_set(self.r, state)
        # cache[self] = result = set(x for x in objects if objects == [y for y in objects if (x, y) not in ext_r or y in ext_c])
        result = set()
        for x in universe:  # TODO COULD BE OPTIMIZED, E.G. IF R HAS EMPTY EXTENSION, ETC.
            ys = ext_c.union(y for y in universe if (x, y) not in ext_r)
            if len(universe) == len(ys):  # No need to compare the sets, as objects has max possible length
                result.add(x)
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Forall(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r.flatten() + self.c.flatten()


class EqualConcept(Concept):
    def __init__(self, r1, r2, sort):
        assert isinstance(r1, Role)
        assert isinstance(r2, Role)
        Concept.__init__(self, sort, 1 + r1.depth + r2.depth)
        self.r1 = r1
        self.r2 = r2
        self.hash = hash((self.__class__, self.r1, self.r2))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.r1 == other.r1 and
                self.r2 == other.r2)

    def extension(self, cache, state, substitution):
        universe = cache.universe_set
        ext_r1 = cache.as_set(self.r1, state)
        ext_r2 = cache.as_set(self.r2, state)

        # cache[self] = result = [x for x in objects if [z for (y, z) in ext_r1 if y == x] == [z for (y, z) in ext_r2 if y == x]]
        result = set()
        for x in universe:
            left = set(z for (y, z) in ext_r1 if y == x)
            right = set(z for (y, z) in ext_r2 if y == x)
            if left == right:
                result.add(x)
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Equal(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r1.flatten() + self.r2.flatten()


class BasicRole(Role):
    def __init__(self, predicate):
        assert isinstance(predicate, tsk.syntax.Predicate)
        super().__init__(predicate.sort, 0)
        self.predicate = predicate

        # This is a bit aggressive, but we assume that predicate names are unique
        self.hash = hash((self.__class__, self.name))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.name == other.name)

    @property
    def name(self):
        return self.predicate.symbol

    def extension(self, cache, state, substitution):
        return cache.as_bitarray(self, state)

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__

    def flatten(self):
        return [self]


class InverseRole(Role):
    def __init__(self, r):
        assert isinstance(r, Role)
        s1, s2 = r.sort
        super().__init__([s2, s1], 1 + r.depth)
        self.r = r
        self.hash = hash((self.__class__, self.r))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.r == other.r)

    def extension(self, cache, state, substitution):
        ext_r = cache.as_set(self.r, state)
        result = set((y, x) for (x, y) in ext_r)
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Inverse(%s)' % repr(self.r)

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r.flatten()


class StarRole(Role):
    def __init__(self, r):
        assert isinstance(r, Role)
        Role.__init__(self, r.sort, 1 + r.depth)
        self.r = r
        self.hash = hash((self.__class__, self.r))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.r == other.r)

    def extension(self, cache, state, substitution):
        ext_r = cache.as_set(self.r, state)
        result = set(transitive_closure(ext_r))
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Star(%s)' % repr(self.r)

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r.flatten()


class CompositionRole(Role):
    def __init__(self, r1, r2):
        assert isinstance(r1, Role)
        assert isinstance(r2, Role)
        Role.__init__(self, [r1.sort[0], r2.sort[1]], 1 + r1.depth + r2.depth)
        self.r1 = r1
        self.r2 = r2
        self.hash = hash((self.__class__, self.r1, self.r2))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.r1 == other.r1 and
                self.r2 == other.r2)

    def extension(self, cache, state, substitution):
        ext_r1 = cache.as_set(self.r1, state)
        ext_r2 = cache.as_set(self.r2, state)
        result = set()
        for a, b in ext_r1:
            for x, y in ext_r2:
                if b == x:
                    result.add((a, y))
                    break  # i.e. break the inner loop
        # for (x, u) in ext_r1:
        #     result.extend((x, z) for (y, z) in ext_r2 if u == y)

        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Composition(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r1.flatten() + self.r2.flatten()


class RestrictRole(Role):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        Role.__init__(self, r.sort, 1 + r.depth + c.depth)
        self.r = r
        self.c = c
        self.hash = hash((self.__class__, self.r, self.c))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.c == other.c and
                self.r == other.r)

    def extension(self, cache, state, substitution):
        ext_c = cache.as_set(self.c, state)
        ext_r = cache.as_set(self.r, state)
        result = set((x, y) for (x, y) in ext_r if y in ext_c)
        return cache.compress(result, self.ARITY)

    def __repr__(self):
        return 'Restrict(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__

    def flatten(self):
        return [self] + self.r.flatten() + self.c.flatten()


class FeatureValueChange(Enum):
    ADD = 1
    DEL = 2
    INC = 3
    DEC = 4
    NIL = 5


class Feature(object):
    def value(self, cache, state, substitution):
        raise NotImplementedError()

    def diff(self, x, y):
        raise NotImplementedError()

    @staticmethod
    def bool_value(value):
        assert value >= 0
        return value > 0

# class NonEmptyConceptFeature(Feature):
#     def __init__(self, c):
#         assert isinstance(c, Concept)
#         super().__init__()
#         self.c = c
#
#     def value(self, cache, state, substitution):
#         """ The feature is true iff the extension of the represented concept has non-empty cardinality"""
#         ext = self.c.extension(cache, state, substitution)
#         return ext.any()
#
#     def diff(self, x, y):
#         assert type(x) is bool
#         assert type(y) is bool
#         if x == y:
#             return FeatureValueChange.NIL
#         if x is True and y is False:
#             return FeatureValueChange.DEL
#         if x is False and y is True:
#             return FeatureValueChange.ADD
#
#     def __repr__(self):
#         return 'NonEmpty({})'.format(self.c)
#
#     __str__ = __repr__


def compute_int_feature_diff(x, y):
    assert type(x) is int and x >= 0
    assert type(y) is int and y >= 0
    if x == y:
        return FeatureValueChange.NIL
    if x > y:
        return FeatureValueChange.DEC
    else:
        return FeatureValueChange.INC


class ConceptCardinalityFeature(Feature):
    def __init__(self, c):
        assert isinstance(c, Concept)
        self.c = c

    def value(self, cache, state, substitution):
        """ The feature value _is_ the cardinality of the extension of the represented concept"""
        ext = self.c.extension(cache, state, substitution)
        return ext.count()

    def diff(self, x, y):
        return compute_int_feature_diff(x, y)

    def __repr__(self):
        return 'cardinality[{}]'.format(self.c)

    __str__ = __repr__


class MinDistanceFeature(Feature):
    def __init__(self, c1, r, c2):
        assert isinstance(c1, Concept) and isinstance(r, Role) and isinstance(c2, Concept)
        self.c1 = c1
        self.r = r
        self.c2 = c2

    def value(self, cache, state, substitution):
        """ The value of the feature is the min distance between any object in the extension of c1 and any object
            on the extension of c2, moving only along r-edges.
        """
        ext_c1 = self.c1.extension(cache, state, substitution)
        ext_c2 = self.c1.extension(cache, state, substitution)
        ext_r = self.r.extension(cache, state, substitution)
        return compute_min_distance(cache.uncompress(ext_c1, self.c1.ARITY),
                                    cache.uncompress(ext_r, self.r.ARITY),
                                    cache.uncompress(ext_c2, self.c2.ARITY),
                                    )

    def diff(self, x, y):
        return compute_int_feature_diff(x, y)

    def __repr__(self):
        return 'min-distance[{}, {}, {}]'.format(self.c1, self.r, self.c2)

    __str__ = __repr__


def most_restricted_type(language, t1, t2):
    if language.is_subtype(t1, t2):
        return t1
    elif language.is_subtype(t2, t1):
        return t2
    return None
