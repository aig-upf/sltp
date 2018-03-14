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
import resource
import argparse
import logging
import sys
import itertools
from signal import signal, SIGPIPE, SIG_DFL

import os
import tarski as tsk
from tarski.io import FstripsReader
from tarski.syntax import builtins

from extensions import UniverseIndex, create_extension_trace
from syntax import Concept, Role, Atom, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, InverseRole, StarRole, RestrictRole, BooleanFeature, Numerical1Feature, Numerical2Feature, \
    AndConcept, most_restricted_type, EmptyConcept
from transitions import read_transitions

signal(SIGPIPE, SIG_DFL)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default=0,
                        type=int)
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    parser.add_argument('--debug', default=False, help='Print additional debugging information')
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
        self.top = UniversalConcept(self.universe_sort)
        self.bot = EmptyConcept(self.universe_sort)
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
        return self.termbox.primitive_concepts, self.termbox.primitive_roles

    def create_atomic_terms(self):
        print('\nConstructing atomic concepts and roles...')

        self.termbox.atomic_concepts = self.create_atomic_concepts(self.termbox.primitive_concepts)
        self.termbox.atomic_roles = self.create_atomic_roles(self.termbox.primitive_roles)

        return self.termbox.atomic_concepts, self.termbox.atomic_roles

    # derive new concepts (1-step derivation) from given set of concepts and roles
    def create_atomic_concepts(self, concepts):
        new_concepts = [self.top, self.bot] + concepts + [self.create_not_concept(c) for c in concepts]
        # TODO Add ParametricConcepts and roles here
        return [term for term in new_concepts if term is not None]

    def create_atomic_roles(self, roles):
        new_roles = list(roles)
        new_roles.extend(InverseRole(r) for r in roles if not isinstance(r, InverseRole))
        new_roles.extend(StarRole(r) for r in roles if not isinstance(r, StarRole))
        new_roles.extend(StarRole(InverseRole(r)) for r in roles if not isinstance(r, InverseRole))
        return [term for term in new_roles if term is not None]

    def derive_concepts(self, old_c, new_c, old_r, new_r):
        result = []

        #
        for fun in (self.create_exists_concept, self.create_forall_concept):
            result.extend(fun(r, c) for r, c in itertools.product(new_r, old_c))
            result.extend(fun(r, c) for r, c in itertools.product(old_r, new_c))
            result.extend(fun(r, c) for r, c in itertools.product(new_r, new_c))

        #
        for fun in (self.create_equal_concept, ):
            result.extend(fun(r1, r2) for r1, r2 in itertools.product(old_r, new_r))
            result.extend(fun(r1, r2) for r1, r2 in itertools.combinations(new_r, 2))

        #
        # for fun in (self.create_and_concept,):
        #     result.extend(fun(c1, c2) for c1, c2 in itertools.product(new_c, old_c))
        #     result.extend(fun(c1, c2) for c1, c2 in itertools.combinations(new_c, 2))

        return [term for term in result if term is not None]

    def derive_roles(self, old_c, new_c, old_r, new_r):
        result = []
        # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
        # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])
        # result.extend([ CompositionRole(r, c) for r, c in itertools.product(roles, roles) if p[0] != p[1] ])

        for fun in (self.create_restrict_role, ):
            result.extend(fun(r, c) for r, c in itertools.product(new_r, old_c))
            result.extend(fun(r, c) for r, c in itertools.product(old_r, new_c))
            result.extend(fun(r, c) for r, c in itertools.product(new_r, new_c))

        return [term for term in result if term is not None]

    def create_exists_concept(self, role: Role, concept: Concept):
        # exists(R.C) = { x | exists y R(x,y) and C(y) }

        result = ExistsConcept(role, concept)
        s1, s2 = role.sort

        if concept == self.bot:
            logging.debug('Concept "{}" is statically empty'.format(result))
            return None

        # TODO ADD: If C is a sort-concept of the same sort than s2, then the concept will be equiv to exist(R.True)
        if not self.language.are_vertically_related(s2, concept.sort):
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
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
            # logging.debug('Concept "{}" equivalent to simpler concept'.format(result))
            return None

        if not self.language.are_vertically_related(s2, concept.sort):
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_and_concept(self, c1: Concept, c2: Concept):
        # C1 AND C2 = { x | C1(x) AND C2(x) }

        result = AndConcept(c1, c2)

        if c1 == c2:
            return None  # No sense in C and C

        if c1 in (self.top, self.bot) or c2 in (self.top, self.bot):
            logging.debug('Concept "{}" pruned, no sense in AND\'ing with top or bot', result)
            return None

        if most_restricted_type(c1.sort.language, c1.sort, c2.sort) is None:
            # i.e. c1 and c2 are disjoint types
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_equal_concept(self, r1: Role, r2: Role):
        assert isinstance(r1, Role) and isinstance(r2, Role)
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
            logging.debug('Concept "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def create_restrict_role(self, r: Role, c: Concept):

        result = RestrictRole(r, c)
        if not self.language.are_vertically_related(r.sort[1], c.sort):
            logging.debug('Role "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        if isinstance(c, UniversalConcept) or c == self.bot:
            logging.debug('Role "{}" pruned; no sense in restricting to top / bot concepts'.format(result))
            return None

        return result


def derive_features(concepts, roles):
    new_features = [BooleanFeature(c) for c in concepts]
    new_features.extend(Numerical1Feature(c) for c in concepts)
    new_features.extend(Numerical2Feature(c1, r, c2) for c1, r, c2 in itertools.product(concepts, roles, concepts))
    return new_features


class InterpretationSet(object):
    def __init__(self, language, state_samples, top, bot):
        self.language = language
        self.state_samples = state_samples
        self.top = top
        self.bot = bot
        self.universe = self.compute_universe(state_samples)
        self.state_extensions = self.initialize_state_extensions(state_samples)
        self.all_traces = dict()  # a dictionary from extension trace to simplest concept / role achieving it

    def generate_extension_trace(self, concept, arity):
        assert isinstance(concept, (Concept, Role))

        trace = []
        for state in self.state_samples:
            cache = self.state_extensions[state]
            ext = concept.extension(cache[self.top], cache, {})
            # assert isinstance(ext, tuple)
            trace.append(ext)

        return create_extension_trace(self.universe, trace, arity)

    def build_cache_for_state(self, state, universe):
        cache = dict()
        cache['atoms'] = []  # Not sure if we'll need this
        for atom in state:
            assert len(atom) in (1, 2, 3)
            name = atom[0]

            if len(atom) == 1:
                cache['atoms'].append(name)
            elif len(atom) == 2:
                cache.setdefault(name, set()).add(universe.index(atom[1]))
            elif len(atom) == 3:
                t = (universe.index(atom[1]), universe.index(atom[2]))
                cache.setdefault(name, set()).add(t)

        cache[self.top] = universe.as_extension()
        cache[self.bot] = set()
        cache["_index_"] = universe
        return cache

    def compute_universe(self, states):
        """ Iterate through all states and collect all possible PDDL objects """
        universe = UniverseIndex()

        for sid, (_, state) in states.items():
            for atom in state:
                assert len(atom) in (1, 2, 3)
                [universe.add(obj) for obj in atom[1:]]

        universe.finish()  # No more objects possible
        return universe

    def initialize_state_extensions(self, states):
        extensions = dict()  # Concept and role extensions indexed by state
        for sid, (_, state) in states.items():
            extensions[sid] = self.build_cache_for_state(state, self.universe)
        return extensions

    def process_term(self, c, arity, name):
        if c is None:
            return None
        trace = self.generate_extension_trace(c, arity)
        try:
            old = self.all_traces[trace]
            # Another concept with same trace exists, so we prune this one
            logging.debug("{} '{}' has equal extension trace to concept '{}'".format(name, c, old))
            return False
        except KeyError:
            self.all_traces[trace] = c
            return True

    def process_concepts(self, elems):
        return (c for c in elems if self.process_term(c, arity=1, name="Concept"))

    def process_roles(self, elems):
        return (c for c in elems if self.process_term(c, arity=2, name="Role"))


def store_terms(concepts, roles, args):
    os.makedirs('terms', exist_ok=True)
    with open(os.path.join('terms', args.transitions_filename + '.concepts'), 'w') as f:
        f.write("\n".join(map(str, concepts)))
    with open(os.path.join('terms', args.transitions_filename + '.roles'), 'w') as f:
        f.write("\n".join(map(str, roles)))


def main(args):
    reader = FstripsReader()
    reader.parse_domain(args.domain)
    problem = reader.problem
    language = problem.language

    print('\nLoading states and transitions...')
    state_samples = read_transitions(args.transitions)

    factory = TerminologicalFactory(language)

    interpretations = InterpretationSet(language, state_samples, factory.top, factory.bot)

    concepts = []
    roles = []

    # Create the "input" terms, i.e. primitive concepts and roles
    factory.create_primitive_terms_from_language()

    # construct primitive concepts and rules: these are input ones plus some other
    ac, ar = factory.create_atomic_terms()

    concepts.extend(interpretations.process_concepts(ac))
    roles.extend(interpretations.process_roles(ar))

    # construct derived concepts and rules obtained with grammar
    c_i, c_j = 0, len(concepts)
    r_i, r_j = 0, len(roles)
    print('\nConstructing derived concepts and roles using {} iteration(s), starting from {} concepts and {} roles'.
          format(args.k, c_j, r_j))

    for i in range(1, args.k+1):

        old_c, new_c = concepts[0:c_i], concepts[c_i:c_j]
        old_r, new_r = roles[0:r_i], roles[r_i:r_j]

        derived_c = factory.derive_concepts(old_c, new_c, old_r, new_r)
        derived_r = factory.derive_roles(old_c, new_c, old_r, new_r)

        print("it. {}: {} concept(s) and {} role(s) generated".format(i, len(derived_c), len(derived_r)))

        concepts.extend(interpretations.process_concepts(derived_c))
        roles.extend(interpretations.process_roles(derived_r))

        c_i, c_j = c_j, len(concepts)
        r_i, r_j = r_j, len(roles)

        print("\t\tof which {} concept(s) and {} role(s) incorporated".format(c_j-c_i, r_j-r_i))

    store_terms(concepts, roles, args)
    print('Final output: %d concept(s) and %d role(s)' % (len(concepts), len(roles)))
    print('Max. memory usage (MB): {}'.format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))


def configure_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO
    filename = os.path.basename(args.transitions)
    args.transitions_filename = '.'.join(filename.split('.')[:-1])
    filename = os.path.join('logs', args.transitions_filename + '.log')
    logging.basicConfig(filename=filename, filemode='w', level=level)


def bootstrap(arguments):
    args = parse_arguments(arguments)
    configure_logging(args)
    return args


if __name__ == "__main__":
    main(bootstrap(sys.argv[1:]))
