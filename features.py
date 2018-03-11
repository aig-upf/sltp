#!/usr/bin/env python

#  Copyright (C) 2018-<date> Blai Bonet
#
#  Permission is hereby granted to distribute this software for
#  non-commercial research purposes, provided that this copyright
#  notice is included with any such distribution.
#
#  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
#  EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
#  SOFTWARE IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
#  ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#  Blai Bonet, bonet@ldc.usb.ve, bonetblai@gmail.com
import sys
import json
import argparse
from signal import signal, SIGPIPE, SIG_DFL
from itertools import product as iter_product
from itertools import combinations as iter_combinations

from tarski.io import FstripsReader
from tarski.syntax import builtins

from utils import check_all_equal, compute_transitive_closure, transitive_closure

signal(SIGPIPE, SIG_DFL)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default='0')
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    return parser.parse_args(args)


# read file line by line
def read_file(filename):
    with open(filename) as f:
        for line in f:
            yield line.rstrip('\n')


# abstract classes for concepts and roles
class Concept(object):
    def __init__(self, type_, depth):
        self.type = type_
        self.depth = depth

    def extension(self, objects, cache, parameter_subst):
        assert False, 'Abstract method should be implemented in derived class'


class Role(object):
    def __init__(self, type_, depth):
        self.type = type_
        self.depth = depth

    def extension(self, objects, cache, parameter_subst):
        assert False, 'Abstract method should be implemented in derived class'


# atoms
class Atom:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


# derived concepts
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
    def __init__(self, name):
        Concept.__init__(self, 'basic', 0)
        self.name = name

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


# derived roles
class BasicRole(Role):
    def __init__(self, name):
        Role.__init__(self, 'basic', 0)
        self.name = name

    def extension(self, objects, cache, parameter_subst):
        return [] if self.name not in cache else cache[self.name]

    def __repr__(self):
        return '%s' % self.name

    __str__ = __repr__


class InverseRole(Role):
    def __init__(self, r):
        assert isinstance(r, Role)
        Role.__init__(self, 'inverse', 1 + r.depth)
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
        Role.__init__(self, 'star', 1 + r.depth)
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
        Role.__init__(self, 'composition', 1 + r1.depth + r2.depth)
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
        Role.__init__(self, 'restrict', 1 + r.depth + c.depth)
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


# derive new concepts (1-step derivation) from given set of concepts and roles

def derive_primitive_concepts(concepts):
    new_concepts = list(concepts)
    new_concepts.extend([NotConcept(c) for c in concepts if not isinstance(c, NotConcept)])
    return new_concepts


def derive_primitive_roles(roles):
    new_roles = list(roles)
    new_roles.extend([InverseRole(r) for r in roles if not isinstance(r, InverseRole)])
    new_roles.extend([StarRole(r) for r in roles if not isinstance(r, StarRole)])
    new_roles.extend([StarRole(InverseRole(r)) for r in roles if not isinstance(r, InverseRole)])
    return new_roles


def derive_concepts(concepts, roles):
    result = []
    # result.extend([ NotConcept(c) for c in concepts if not isinstance(c, NotConcept) and c.depth > 0 ])
    # result.extend([ NotConcept(c) for c in concepts if not isinstance(c, NotConcept) ])
    # result.extend(map(lambda p: AndConcept(p[0], p[1]), iter_combinations(concepts, 2)))
    result.extend([ExistsConcept(p[0], p[1]) for p in iter_product(roles, concepts)])
    result.extend([ForallConcept(p[0], p[1]) for p in iter_product(roles, concepts)])
    result.extend([EqualConcept(p[0], p[1]) for p in iter_combinations(roles, 2)])
    return result


def derive_roles(concepts, roles):
    result = []
    # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
    # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])
    # result.extend([ CompositionRole(p[0], p[1]) for p in iter_product(roles, roles) if p[0] != p[1] ])
    result.extend([RestrictRole(p[0], p[1]) for p in iter_product(roles, concepts)])
    return result


def extend_concepts_and_roles(concepts, roles):
    new_concepts = derive_concepts(concepts, roles)
    new_roles = derive_roles(concepts, roles)
    return new_concepts, new_roles


# features
class Feature(object):
    def __init__(self, type_):
        self.type = type_


class BooleanFeature(Feature):
    def __init__(self, c):
        assert isinstance(c, Concept)
        Feature.__init__(self, 'boolean')
        self.c = c

    def __repr__(self):
        return 'Boolean(%s)' % repr(self.c)

    __str__ = __repr__


class Numerical1Feature(Feature):
    def __init__(self, c):
        assert isinstance(c, Concept)
        Feature.__init__(self, 'boolean')
        self.c = c

    def __repr__(self):
        return 'Numerical1(%s)' % repr(self.c)

    __str__ = __repr__


class Numerical2Feature(Feature):
    def __init__(self, c1, r, c2):
        assert isinstance(c1, Concept)
        assert isinstance(r, Role)
        assert isinstance(c2, Concept)
        Feature.__init__(self, 'boolean')
        self.c1 = c1
        self.r = r
        self.c2 = c2

    def __repr__(self):
        return 'Numerical2(%s,%s,%s)' % repr(self.c1, self.r, self.c2)

    __str__ = __repr__


def build_cache_for_state(state):
    cache = {'<universe>': set(), 'atoms': []}
    for atom in state:
        assert atom
        name = atom[0]

        if len(atom) == 1:
            cache['atoms'].append(name)
        else:
            if name not in cache:
                cache[name] = []
            if len(atom) == 2:
                cache[name].append(atom[1])
                cache['<universe>'].add(atom[1])
            elif len(atom) == 3:
                cache[name].append((atom[1], atom[2]))
                cache['<universe>'].add(atom[1])
                cache['<universe>'].add(atom[2])
            else:
                assert False, "Expecting at most two arguments, but got '(%s)'" % ' '.join(atom)
    return cache


def derive_features(concepts, roles):
    new_features = [BooleanFeature(c) for c in concepts]
    new_features.extend(Numerical1Feature(c) for c in concepts)
    new_features.extend(Numerical2Feature(p[0], p[1], p[2]) for p in iter_product(concepts, roles, concepts))
    return new_features


# read signature
g_map_signature = {}
g_input_atoms = []
g_input_concepts = []
g_input_roles = []


def read_signature(language):
    assert len(language.functions) == 0

    for predicate in language.predicates:
        if builtins.is_builtin_predicate(predicate):
            continue  # Skip "=" and other built-in symbols

        name = predicate.symbol
        g_map_signature[name] = predicate.arity

        if predicate.arity == 0:
            g_input_atoms.append(Atom(name))
        elif predicate.arity == 1:
            g_input_concepts.append(BasicConcept(name))
        elif predicate.arity == 2:
            g_input_roles.append(BasicRole(name))
        else:
            print("Predicate {} with arity > 2 ignored".format(predicate))

    print('%d input atoms (0-ary signature):' % len(g_input_atoms), [str(item) for item in g_input_atoms])
    print('%d input concepts (1-ary signature):' % len(g_input_concepts), [str(item) for item in g_input_concepts])
    print('%d input roles (2-ary signature):' % len(g_input_roles), [str(item) for item in g_input_roles])


# read transitions
g_states_by_str = {}
g_states_by_id = {}
g_transitions = {}


def read_transitions(transitions_filename):
    raw_file = [line.replace(' ', '') for line in read_file(transitions_filename) if line[0:6] == '{"id":']
    for raw_line in raw_file:
        j = json.loads(raw_line)
        j_atoms = [atom.rstrip(')').replace('(', ',').split(',') for atom in j['atoms']]
        j_atoms_str = str(j_atoms)
        j_id = -1

        # insert state into hash with (normalized) id
        if j_atoms_str not in g_states_by_str:
            j_id = int(j['id'])
            g_states_by_str[j_atoms_str] = g_states_by_id[j_id] = (j_id, j_atoms)
        else:
            j_id = g_states_by_str[j_atoms_str][0]

        # insert (normalized) transition into hash
        if j['parent'] != j['id']:
            j_pid = int(j['parent'])
            jp = json.loads(raw_file[j_pid])
            assert jp['id'] == j_pid
            jp_atoms = [atom.rstrip(')').replace('(', ',').split(',') for atom in jp['atoms']]
            jp_atoms_str = str(jp_atoms)

            assert jp_atoms_str in g_states_by_str
            jp_id = g_states_by_str[jp_atoms_str][0]
            if jp_id not in g_transitions: g_transitions[jp_id] = []
            g_transitions[jp_id].append(j_id)

    # check soundness
    for src in g_transitions:
        assert src in g_states_by_id
        for dst in g_transitions[src]:
            assert dst in g_states_by_id

    print('#lines-raw-file=%d, #state-by-str=%d, #states-by-id=%d, #transition-entries=%d, #transitions=%d' % (
        len(raw_file), len(g_states_by_str), len(g_states_by_id), len(g_transitions),
        sum([len(g_transitions[src]) for src in g_transitions])))


def prune_elements(existing, novel, extensions):

    assert extensions

    good = []
    for element in novel:
        all_elem_extensions = [element.extension(cache['<universe>'], cache, {})
                               for _, cache in extensions.items()]

        # if check_all_equal(all_elem_extensions):
        #     print('Concept/role "{}" has constant extension over all samples'.format(element))

        if check_all_equal(all_elem_extensions) and len(all_elem_extensions[1]) == 0:
            # print('Concept/role "{}" has empty extension over all samples'.format(element))
            pass
        else:
            good.append(element)

        # for state_id, state_cache in extensions.items():
        #     elem_ext = element.extension(state_cache['<universe>'], state_cache, {})

    print("{} concepts/roles pruned because of emptyness".format(len(novel) - len(good)))

    # l = list(itertools.product(existing, novel))
    return good


def print_atstar(state_extensions):
    s = state_extensions[0]
    at_closure = StarRole(BasicRole('at'))
    e = at_closure.extension(s['<universe>'], s, {})
    print("Extension of AT*: {}".format(e))


def main(args):
    reader = FstripsReader()
    reader.parse_domain(args.domain)
    problem = reader.problem

    print('Loading domain predicate signatures...')
    read_signature(problem.language)

    print('\nLoading states and transitions...')
    read_transitions(args.transitions)

    # construct primitive concepts and rules: these are input ones plus some other
    print('\nConstructing primitive concepts and roles...')
    g_primitive_concepts = derive_primitive_concepts(g_input_concepts)
    g_primitive_roles = derive_primitive_roles(g_input_roles)
    print('%d primitive concept(s):' % len(g_primitive_concepts),
          [(str(item), item.depth) for item in g_primitive_concepts])
    print('%d primitive rule(s):' % len(g_primitive_roles), [(str(item), item.depth) for item in g_primitive_roles])

    state_extensions = dict()  # Concept and role extensions indexed by state
    update_state_extensions(g_primitive_concepts + g_primitive_roles, state_extensions)
    g_primitive_concepts = prune_elements([], g_primitive_concepts, state_extensions)
    g_primitive_roles = prune_elements([], g_primitive_roles, state_extensions)

    # construct derived concepts and rules obtained with grammar
    print('\nConstructing derived concepts and roles using %d iteration(s)...' % int(args.k))
    print_atstar(state_extensions)
    g_derived_concepts = list(g_primitive_concepts)
    g_derived_roles = list(g_primitive_roles)
    for i in range(0, int(args.k)):
        new_concepts, new_roles = extend_concepts_and_roles(g_derived_concepts, g_derived_roles)
        print('iteration %d: %d new concept(s) and %d new role(s)' % (1 + i, len(new_concepts), len(new_roles)))

        update_state_extensions(new_concepts + new_roles, state_extensions)
        new_concepts = prune_elements(g_derived_concepts, new_concepts, state_extensions)
        new_roles = prune_elements(g_derived_roles, new_roles, state_extensions)

        g_derived_concepts.extend(new_concepts)
        g_derived_roles.extend(new_roles)

    print('final: %d concept(s) and %d role(s)' % (len(g_derived_concepts), len(g_derived_roles)))
    print('\n\nDerived concepts:\n' + '\n'.join("{} ({})".format(item, item.depth) for item in g_derived_concepts))
    print('\n\nDerived roles:\n' + '\n'.join("{} ({})".format(item, item.depth) for item in g_derived_roles))

    # Tests
    print_atstar(state_extensions)



def update_state_extensions(elements, state_extensions):
    """ Updates the given state extensions with the extensions in each state of all
    given elements (concepts and roles)"""
    for sid in g_states_by_id:
        compute_state_extensions(sid, elements, state_extensions)


def compute_state_extensions(state, elements, state_extensions):
    # Retrieve the state cache if previously built, otherwise create a new one
    cache = state_extensions.get(state, None) or build_cache_for_state(g_states_by_id[state][1])
    state_extensions[state] = cache

    for elem in elements:
        assert isinstance(elem, (Concept, Role))
        _ = elem.extension(cache['<universe>'], cache, {})
        # if ext: print('Extension of CONCEPT %s: %s' % (str(c), _))


if __name__ == "__main__":
    main(parse_arguments(sys.argv[1:]))
