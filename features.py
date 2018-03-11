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
import tarski as tsk


from syntax import Concept, Role, Atom, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, InverseRole, StarRole, RestrictRole, BooleanFeature, Numerical1Feature, Numerical2Feature
from utils import check_all_equal, read_file

signal(SIGPIPE, SIG_DFL)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default='0')
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    return parser.parse_args(args)


class TermBox(object):
    """ A simple holder object where we'll classify the different concepts and roles that we generate """

    def __init__(self):
        self.primitive_atoms = []
        self.primitive_concepts = []
        self.primitive_roles = []

        self.atomic_concepts = []
        self.atomic_roles = []

    def all_primitives(self):
        return self.primitive_concepts + self.primitive_roles

    def all_atoms(self):
        return self.atomic_concepts + self.atomic_roles

    def print_primitives(self):
        print('%d input atoms (0-ary signature):' % len(self.primitive_atoms), [str(item) for item in self.primitive_atoms])
        print('%d input concepts (1-ary signature):' % len(self.primitive_concepts), [str(item) for item in self.primitive_concepts])
        print('%d input roles (2-ary signature):' % len(self.primitive_roles), [str(item) for item in self.primitive_roles])

    def print_atoms(self):
        print('%d primitive concept(s):' % len(self.atomic_concepts),
              [(str(item), item.depth) for item in self.atomic_concepts])
        print('%d primitive rule(s):' % len(self.atomic_roles), [(str(item), item.depth) for item in self.atomic_roles])


class TerminologicalFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage):
        self.language = language
        self.universe_sort = language.get_sort('object')
        self.termbox = TermBox()

        self.map_signature = {}

    def create_primitive_terms_from_language(self):
        assert len(self.language.functions) == 0
        for predicate in self.language.predicates:
            if builtins.is_builtin_predicate(predicate):
                continue  # Skip "=" and other built-in symbols

            name = predicate.symbol
            self.map_signature[name] = predicate.arity

            if predicate.arity == 0:
                self.termbox.primitive_atoms.append(Atom(name))
            elif predicate.arity == 1:
                self.termbox.primitive_concepts.append(BasicConcept(predicate))
            elif predicate.arity == 2:
                self.termbox.primitive_roles.append(BasicRole(predicate))
            else:
                print("Predicate {} with arity > 2 ignored".format(predicate))

        self.termbox.print_primitives()

    def create_atomic_terms(self):
        print('\nConstructing atomic concepts and roles...')

        self.termbox.atomic_concepts = self.derive_atomic_concepts(self.termbox.primitive_concepts)
        self.termbox.atomic_roles = self.derive_atomic_roles(self.termbox.primitive_roles)

        # g_primitive_concepts = self.derive_atomic_concepts(self.termbox.primitive_concepts)
        # g_primitive_roles = self.derive_atomic_roles(self.termbox.primitive_roles)

    # derive new concepts (1-step derivation) from given set of concepts and roles
    def derive_atomic_concepts(self, concepts):
        new_concepts = list(concepts)
        universal = UniversalConcept(self.universe_sort)
        new_concepts.append(universal)
        new_concepts.append(self.create_not_concept(universal))
        new_concepts.extend(self.create_not_concept(c) for c in concepts)
        # TODO Add ParametricConcepts and roles here
        return [term for term in new_concepts if term is not None]

    def derive_atomic_roles(self, roles):
        new_roles = list(roles)
        new_roles.extend(InverseRole(r) for r in roles if not isinstance(r, InverseRole))
        new_roles.extend(StarRole(r) for r in roles if not isinstance(r, StarRole))
        new_roles.extend(StarRole(InverseRole(r)) for r in roles if not isinstance(r, InverseRole))
        return [term for term in new_roles if term is not None]


    def derive_concepts(self, concepts, roles):
        result = []
        # result.extend([ NotConcept(c) for c in concepts if not isinstance(c, NotConcept) and c.depth > 0 ])
        # result.extend([ NotConcept(c) for c in concepts if not isinstance(c, NotConcept) ])
        # result.extend(map(lambda p: AndConcept(r, c), iter_combinations(concepts, 2)))
        result.extend(self.create_exists_concept(r, c) for r, c in iter_product(roles, concepts))
        result.extend(self.create_forall_concept(r, c) for r, c in iter_product(roles, concepts))
        result.extend(self.create_equal_concept(r1, r2) for r1, r2 in iter_combinations(roles, 2))
        return [term for term in result if term is not None]


    def derive_roles(self, concepts, roles):
        result = []
        # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
        # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])
        # result.extend([ CompositionRole(r, c) for r, c in iter_product(roles, roles) if p[0] != p[1] ])
        result.extend(RestrictRole(r, c) for r, c in iter_product(roles, concepts))
        return [term for term in result if term is not None]


    def extend_concepts_and_roles(self, concepts, roles):
        new_concepts = self.derive_concepts(concepts, roles)
        new_roles = self.derive_roles(concepts, roles)
        return new_concepts, new_roles

    def create_exists_concept(self, role: Role, concept: Concept):
        # exists(R.C) = { x | exists y R(x,y) and C(y) }

        result = ExistsConcept(role, concept)
        s1, s2 = role.sort

        # TODO ADD: If C is a sort-concept of the same sort than s2, then the concept will be equiv to exist(R.True)
        if not self.language.are_vertically_related(s2, concept.sort):
            print('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_not_concept(self, concept: Concept):
        if isinstance(concept, NotConcept):
            return None

        result = NotConcept(concept, self.universe_sort)
        return result

    def create_forall_concept(self, role: Role, concept: Concept):
        # forall(R.C) = { x | forall y R(x,y) implies C(y) }

        result = ForallConcept(role, concept)
        s1, s2 = role.sort

        if isinstance(concept, UniversalConcept):
            # print('Concept "{}" equivalent to simpler concept'.format(result))
            return None

        if not self.language.are_vertically_related(s2, concept.sort):
            print('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_equal_concept(self, r1: Role, r2: Role):

        # The sort of the resulting concept will be the most restricted sort between the left sorts of the two roles
        if self.language.is_subtype(r1.sort[0], r2.sort[0]):
            sort = r1.sort[0]
        elif self.language.is_subtype(r2.sort[0], r1.sort[0]):
            sort = r2.sort[0]
        else:
            sort = None

        result = EqualConcept(r1, r2, sort)

        if not self.language.are_vertically_related(r1.sort[0], r2.sort[0]) or \
                not self.language.are_vertically_related(r1.sort[1], r2.sort[1]):
            print('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result


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


# read transitions
g_states_by_str = {}
g_states_by_id = {}
g_transitions = {}


def normalize_atom_name(name):
    return name.replace('()', '').rstrip(')').replace('(', ',').split(',')


def read_transitions(transitions_filename):
    raw_file = [line.replace(' ', '') for line in read_file(transitions_filename) if line[0:6] == '{"id":']
    for raw_line in raw_file:
        j = json.loads(raw_line)
        j_atoms = [normalize_atom_name(atom) for atom in j['atoms']]
        j_atoms_str = str(j_atoms)

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
            jp_atoms = [normalize_atom_name(atom) for atom in jp['atoms']]

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


def prune_elements(existing, novel, extensions, name):

    assert extensions

    all_universes = [cache['<universe>'] for cache in extensions.values()]
    assert check_all_equal(all_universes), "We should expect the universe of interpretation in all states to be equal"
    universe = all_universes[0]

    good = []
    for element in novel:
        all_elem_extensions = [element.extension(cache['<universe>'], cache, {})
                               for cache in extensions.values()]

        # if check_all_equal(all_elem_extensions):
        #     print('Concept/role "{}" has constant extension over all samples'.format(element))

        all_equal = check_all_equal(all_elem_extensions)

        if all_equal:
            if repr(element) != "Not(<universe>)" and len(all_elem_extensions[1]) == 0:
                print('{} "{}" has empty extension over all samples'.format(name, element))
                continue

            if repr(element) != "<universe>" and len(all_elem_extensions[1]) == len(universe):
                print('{} "{}" has universal extension over all samples'.format(name, element))
                continue

        good.append(element)

        # for state_id, state_cache in extensions.items():
        #     elem_ext = element.extension(state_cache['<universe>'], state_cache, {})

    print("{} {}(s) pruned because of emptyness".format(len(novel) - len(good), name))

    # l = list(itertools.product(existing, novel))
    return good


def main(args):
    reader = FstripsReader()
    reader.parse_domain(args.domain)
    problem = reader.problem
    language = problem.language

    print('\nLoading states and transitions...')
    read_transitions(args.transitions)

    factory = TerminologicalFactory(language)
    termbox = factory.termbox

    # Create the "input" terms, i.e. primitive concepts and roles
    factory.create_primitive_terms_from_language()

    # construct primitive concepts and rules: these are input ones plus some other
    factory.create_atomic_terms()

    state_extensions = dict()  # Concept and role extensions indexed by state
    update_state_extensions(termbox.all_atoms(), state_extensions)
    termbox.atomic_concepts = prune_elements([], termbox.atomic_concepts, state_extensions, "concept")
    termbox.atomic_roles = prune_elements([], termbox.atomic_roles, state_extensions, "role")

    # construct derived concepts and rules obtained with grammar
    print('\nConstructing derived concepts and roles using %d iteration(s)...' % int(args.k))
    derived_concepts = list(termbox.atomic_concepts)
    derived_roles = list(termbox.atomic_roles)
    for i in range(0, int(args.k)):
        new_concepts, new_roles = factory.extend_concepts_and_roles(derived_concepts, derived_roles)
        print('iteration %d: %d new concept(s) and %d new role(s)' % (1 + i, len(new_concepts), len(new_roles)))

        update_state_extensions(new_concepts + new_roles, state_extensions)
        new_concepts = prune_elements(derived_concepts, new_concepts, state_extensions, "concept")
        new_roles = prune_elements(derived_roles, new_roles, state_extensions, "role")

        derived_concepts.extend(new_concepts)
        derived_roles.extend(new_roles)

    print('final: %d concept(s) and %d role(s)' % (len(derived_concepts), len(derived_roles)))
    # print('\n\nDerived concepts:\n' + '\n'.join("{} ({})".format(item, item.depth) for item in derived_concepts))
    # print('\n\nDerived roles:\n' + '\n'.join("{} ({})".format(item, item.depth) for item in derived_roles))


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
