import tarski as tsk

from utils import transitive_closure


# abstract classes for concepts and roles
class Concept(object):
    def __init__(self, sort, depth):
        self.sort = sort
        self.depth = depth

    def extension(self, objects, cache, parameter_subst):
        assert False, 'Abstract method should be implemented in derived class'


class Role(object):
    def __init__(self, sort, depth):
        self.sort = sort
        self.depth = depth

    def extension(self, objects, cache, parameter_subst):
        assert False, 'Abstract method should be implemented in derived class'


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

    def extension(self, objects, cache, parameter_subst):
        return cache[self]

    def __repr__(self):
        return '<universe>'

    __str__ = __repr__


class EmptyConcept(Concept):
    def __init__(self, universal_sort):
        Concept.__init__(self, universal_sort, 0)
        self.hash = hash(self.__class__)

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return self.__class__ is other.__class__

    def extension(self, objects, cache, parameter_subst):
        return set()

    def __repr__(self):
        return '<empty>'

    __str__ = __repr__


class ParametricConcept(Concept):
    def __init__(self, parameter):
        Concept.__init__(self, 'parametric', 0)
        self.parameter = parameter
        self.hash = hash((self.__class__, self.parameter))

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        return (hasattr(other, 'hash') and self.hash == other.hash and self.__class__ is other.__class__ and
                self.parameter == other.parameter)

    def extension(self, objects, cache, parameter_subst):
        assert self.parameter in parameter_subst, "Parameter '?%s' should appear in substitution: %s" % (
            self.parameter, str(parameter_subst))
        assert False, "Revise implementation details of following line!"
        return set(parameter_subst[self.parameter])

    def __repr__(self):
        return '?%s' % self.parameter

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        return set() if self.name not in cache else cache[self.name]

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_c = self.c.extension(objects, cache, parameter_subst)
            cache[self] = result = objects - ext_c
            return result

    def __repr__(self):
        return 'Not(%s)' % repr(self.c)

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_c1 = self.c1.extension(objects, cache, parameter_subst)
            ext_c2 = self.c2.extension(objects, cache, parameter_subst)
            cache[self] = result = ext_c1 & ext_c2
            return result

    def __repr__(self):
        return 'And(%s,%s)' % (repr(self.c1), repr(self.c2))

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            # result = [x for x in objects if [z for (y, z) in ext_r if y == x and z in ext_c]]
            cache[self] = result = set(x for x, y in ext_r if y in ext_c)
            return result

    def __repr__(self):
        return 'Exists(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r = self.c.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            # cache[self] = result = set(x for x in objects if objects == [y for y in objects if (x, y) not in ext_r or y in ext_c])
            cache[self] = result = set()
            for x in objects:
                ys = ext_c.union(y for y in objects if (x, y) not in ext_r)
                if len(objects) == len(ys):  # No need to compare the sets, as objects has max possible length
                    result.add(x)
            return result

    def __repr__(self):
        return 'Forall(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r1 = self.r1.extension(objects, cache, parameter_subst)
            ext_r2 = self.r2.extension(objects, cache, parameter_subst)
            # cache[self] = result = [x for x in objects if [z for (y, z) in ext_r1 if y == x] == [z for (y, z) in ext_r2 if y == x]]
            cache[self] = result = set()
            for x in objects:
                left = set(z for (y, z) in ext_r1 if y == x)
                right = set(z for (y, z) in ext_r2 if y == x)
                if left == right:
                    result.add(x)
            return result

    def __repr__(self):
        return 'Equal(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        return set() if self.name not in cache else cache[self.name]

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            cache[self] = result = set((y, x) for (x, y) in ext_r)
            return result

    def __repr__(self):
        return 'Inverse(%s)' % repr(self.r)

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self not in cache:
            ext = self.r.extension(objects, cache, parameter_subst)
            # cache[self] = compute_transitive_closure(ext)
            cache[self] = result = set(transitive_closure(ext))

        return cache[self]

    def __repr__(self):
        return 'Star(%s)' % repr(self.r)

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r1 = self.r1.extension(objects, cache, parameter_subst)
            ext_r2 = self.r2.extension(objects, cache, parameter_subst)
            cache[self] = result = set()
            for a, b in ext_r1:
                for x, y in ext_r2:
                    if b == x:
                        result.add((a, y))
                        break  # i.e. break the inner loop
            # for (x, u) in ext_r1:
            #     result.extend((x, z) for (y, z) in ext_r2 if u == y)

            return result

    def __repr__(self):
        return 'Composition(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__


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

    def extension(self, objects, cache, parameter_subst):
        if self in cache:
            return cache[self]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            cache[self] = result = set((x, y) for (x, y) in ext_r if y in ext_c)
            return result

    def __repr__(self):
        return 'Restrict(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__


class Feature(object):
    def __init__(self):
        pass


class BooleanFeature(Feature):
    def __init__(self, c):
        assert isinstance(c, Concept)
        super().__init__()
        self.c = c

    def __repr__(self):
        return 'Boolean(%s)' % repr(self.c)

    __str__ = __repr__


class Numerical1Feature(Feature):
    def __init__(self, c):
        assert isinstance(c, Concept)
        super().__init__()
        self.c = c

    def __repr__(self):
        return 'Numerical1(%s)' % repr(self.c)

    __str__ = __repr__


class Numerical2Feature(Feature):
    def __init__(self, c1, r, c2):
        assert isinstance(c1, Concept)
        assert isinstance(r, Role)
        assert isinstance(c2, Concept)
        super().__init__()
        self.c1 = c1
        self.r = r
        self.c2 = c2

    def __repr__(self):
        return 'Numerical2(%s,%s,%s)' % repr(self.c1, self.r, self.c2)

    __str__ = __repr__


def most_restricted_type(language, t1, t2):
    if language.is_subtype(t1, t2):
        return t1
    elif language.is_subtype(t2, t1):
        return t2
    return None
