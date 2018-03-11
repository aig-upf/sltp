#!/usr/bin/env python

import argparse
import json
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
from itertools import combinations as iter_combinations
from itertools import product as iter_product
from signal import signal, SIGPIPE, SIG_DFL

from tarski.io import FstripsReader
from tarski.syntax import builtins

from syntax import Concept, Role, Atom, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, InverseRole, StarRole, RestrictRole
from utils import check_all_equal, read_file

signal(SIGPIPE, SIG_DFL)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default='0')
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    return parser.parse_args(args)


# derive new concepts (1-step derivation) from given set of concepts and roles

def derive_primitive_concepts(concepts):
    new_concepts = list(concepts)
    new_concepts.append(UniversalConcept())
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
    # result.extend(map(lambda p: AndConcept(r, c), iter_combinations(concepts, 2)))
    result.extend([ExistsConcept(r, c) for r, c in iter_product(roles, concepts)])
    result.extend([ForallConcept(r, c) for r, c in iter_product(roles, concepts)])
    result.extend([EqualConcept(r, c) for r, c in iter_combinations(roles, 2)])
    return result


def derive_roles(concepts, roles):
    result = []
    # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
    # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])
    # result.extend([ CompositionRole(r, c) for r, c in iter_product(roles, roles) if p[0] != p[1] ])
    result.extend([RestrictRole(r, c) for r, c in iter_product(roles, concepts)])
    return result


def extend_concepts_and_roles(concepts, roles):
    new_concepts = derive_concepts(concepts, roles)
    new_roles = derive_roles(concepts, roles)
    return new_concepts, new_roles


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
    new_features.extend(Numerical2Feature(c1, r, c2) for c1, r, c2 in iter_product(concepts, roles, concepts))
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
            g_input_concepts.append(BasicConcept(predicate))
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
    language = problem.language

    print('Loading domain predicate signatures...')
    read_signature(language)

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
