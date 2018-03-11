import tarski as tsk

from utils import transitive_closure


# abstract classes for concepts and roles
class Concept(object):
    def __init__(self, sort, depth):
        self.depth = depth
        self.sort = sort

    def extension(self, objects, cache, parameter_subst):
        assert False, 'Abstract method should be implemented in derived class'


class Role(object):
    def __init__(self, depth):
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
    def __init__(self):
        Concept.__init__(self, 'universal', 0)

    def extension(self, objects, cache, parameter_subst):
        return [] if repr(self) not in cache else cache[repr(self)]

    def __repr__(self):
        return '<universe>'

    __str__ = __repr__


class ParametricConcept(Concept):
    def __init__(self, parameter):
        Concept.__init__(self, 'parametric', 0)
        self.parameter = parameter

    def extension(self, objects, cache, parameter_subst):
        assert self.parameter in parameter_subst, "Parameter '?%s' should appear in substitution: %s" % (
            self.parameter, str(parameter_subst))
        return list(parameter_subst[self.parameter])

    def __repr__(self):
        return '?%s' % self.parameter

    __str__ = __repr__


class BasicConcept(Concept):
    def __init__(self, predicate):
        assert isinstance(predicate, tsk.syntax.Predicate)
        Concept.__init__(self, predicate.sort, 0)
        self.predicate = predicate

    @property
    def name(self):
        return self.predicate.symbol

    def extension(self, objects, cache, parameter_subst):
        return [] if self.name not in cache else cache[self.name]

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


class NotConcept(Concept):
    def __init__(self, c):
        assert isinstance(c, Concept)
        Concept.__init__(self, 'not', 1 + c.depth)
        self.c = c

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_c = self.c.extension(objects, cache, parameter_subst)
            result = [x for x in objects if x not in ext_c]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Not(%s)' % repr(self.c)

    __str__ = __repr__


class AndConcept(Concept):
    def __init__(self, c1, c2):
        assert isinstance(c1, Concept)
        assert isinstance(c2, Concept)
        Concept.__init__(self, 'and', 1 + c1.depth + c2.depth)
        self.c1 = c1
        self.c2 = c2

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_c1 = self.c1.extension(objects, cache, parameter_subst)
            ext_c2 = self.c2.extension(objects, cache, parameter_subst)
            result = [x for x in ext_c1 if x in ext_c2]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'And(%s,%s)' % (repr(self.c1), repr(self.c2))

    __str__ = __repr__


class ExistsConcept(Concept):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        Concept.__init__(self, 'exists', 1 + r.depth + c.depth)
        self.r = r
        self.c = c

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            result = [x for x in objects if [z for (y, z) in ext_r if y == x and z in ext_c]]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Exists(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__


class ForallConcept(Concept):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        Concept.__init__(self, 'forall', 1 + r.depth + c.depth)
        self.r = r
        self.c = c

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r = self.c.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            result = [x for x in objects if objects == [y for y in objects if (x, y) not in ext_r or y in ext_c]]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Forall(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__


class EqualConcept(Concept):
    def __init__(self, r1, r2):
        assert isinstance(r1, Role)
        assert isinstance(r2, Role)
        Concept.__init__(self, 'equal', 1 + r1.depth + r2.depth)
        self.r1 = r1
        self.r2 = r2

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r1 = self.r1.extension(objects, cache, parameter_subst)
            ext_r2 = self.r2.extension(objects, cache, parameter_subst)
            result = [x for x in objects if [z for (y, z) in ext_r1 if y == x] == [z for (y, z) in ext_r2 if y == x]]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Equal(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__


class BasicRole(Role):
    def __init__(self, name):
        Role.__init__(self, 0)
        self.name = name

    def extension(self, objects, cache, parameter_subst):
        return [] if self.name not in cache else cache[self.name]

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


class InverseRole(Role):
    def __init__(self, r):
        assert isinstance(r, Role)
        Role.__init__(self, 1 + r.depth)
        self.r = r

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            result = [(y, x) for (x, y) in ext_r]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Inverse(%s)' % repr(self.r)

    __str__ = __repr__


class StarRole(Role):
    def __init__(self, r):
        assert isinstance(r, Role)
        Role.__init__(self, 1 + r.depth)
        self.r = r

    def extension(self, objects, cache, parameter_subst):
        if repr(self) not in cache:
            ext = self.r.extension(objects, cache, parameter_subst)
            # cache[repr(self)] = compute_transitive_closure(ext)
            cache[repr(self)] = transitive_closure(ext)

        return cache[repr(self)]

    def __repr__(self):
        return 'Star(%s)' % repr(self.r)

    __str__ = __repr__


class CompositionRole(Role):
    def __init__(self, r1, r2):
        assert isinstance(r1, Role)
        assert isinstance(r2, Role)
        Role.__init__(self, 1 + r1.depth + r2.depth)
        self.r1 = r1
        self.r2 = r2

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r1 = self.r1.extension(objects, cache, parameter_subst)
            ext_r2 = self.r2.extension(objects, cache, parameter_subst)
            result = []
            for (x, u) in ext_r1:
                result.extend([(x, z) for (y, z) in ext_r2 if u == y])
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Composition(%s,%s)' % (repr(self.r1), repr(self.r2))

    __str__ = __repr__


class RestrictRole(Role):
    def __init__(self, r, c):
        assert isinstance(r, Role)
        assert isinstance(c, Concept)
        Role.__init__(self, 1 + r.depth + c.depth)
        self.r = r
        self.c = c

    def extension(self, objects, cache, parameter_subst):
        if repr(self) in cache:
            return cache[repr(self)]
        else:
            ext_r = self.r.extension(objects, cache, parameter_subst)
            ext_c = self.c.extension(objects, cache, parameter_subst)
            result = [(x, y) for (x, y) in ext_r if y in ext_c]
            cache[repr(self)] = result
            return result

    def __repr__(self):
        return 'Restrict(%s,%s)' % (repr(self.r), repr(self.c))

    __str__ = __repr__



# features
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
