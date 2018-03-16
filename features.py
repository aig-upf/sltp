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
import time

from bitarray import bitarray
from tarski.io import FstripsReader
from tarski.syntax import builtins

# import profiling
from extensions import UniverseIndex, ExtensionCache
from syntax import Concept, Role, Atom, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, InverseRole, StarRole, RestrictRole, BooleanFeature, Numerical1Feature, Numerical2Feature, \
    AndConcept, most_restricted_type, EmptyConcept, CompositionRole
from transitions import read_transitions
signal(SIGPIPE, SIG_DFL)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default=0,
                        type=int)
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args(args)


class TermBox(object):
    """ A simple holder object where we'll classify the different concepts and roles that we generate """

    def __init__(self):
        self.primitive_atoms = []


class TerminologicalFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage):
        self.language = language
        self.universe_sort = language.get_sort('object')
        self.top = UniversalConcept(self.universe_sort)
        self.bot = EmptyConcept(self.universe_sort)
        self.termbox = TermBox()

        self.map_signature = {}

    def create_primitives(self):
        assert len(self.language.functions) == 0

        concepts, roles = [], []
        for predicate in self.language.predicates:
            if builtins.is_builtin_predicate(predicate):
                continue  # Skip "=" and other built-in symbols

            name = predicate.symbol
            self.map_signature[name] = predicate.arity

            if predicate.arity == 0:
                self.termbox.primitive_atoms.append(Atom(name))
            elif predicate.arity == 1:
                concepts.append(BasicConcept(predicate))
            elif predicate.arity == 2:
                roles.append(BasicRole(predicate))
            else:
                print("Predicate {} with arity > 2 ignored".format(predicate))

        atoms = self.termbox.primitive_atoms
        print('{} input (nullary) atoms : {}'.format(len(atoms), list(map(str, atoms))))
        print('{} input (unary) concepts: {}'.format(len(concepts), list(map(str, concepts))))
        print('{} input (binary) roles  : {}'.format(len(roles), list(map(str, roles))))
        return concepts, roles

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
        for fun in (self.create_and_concept,):
            result.extend(fun(c1, c2) for c1, c2 in itertools.product(new_c, old_c))
            result.extend(fun(c1, c2) for c1, c2 in itertools.combinations(new_c, 2))

        return [term for term in result if term is not None]

    def derive_roles(self, old_c, new_c, old_r, new_r):
        result = []
        # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
        # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])

        for fun in (self.create_composition_role, ):
            result.extend(fun(r1, r2) for r1, r2 in itertools.product(old_r, new_r))
            result.extend(fun(r1, r2) for r1, r2 in itertools.combinations(new_r, 2))

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

        if isinstance(role, RestrictRole) and concept == self.top:
            return None  # Is syntactically equivalent to a simpler exists concept

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
            logging.debug('Concept "%s" pruned, no sense in AND\'ing with top or bot', result)
            return None

        if most_restricted_type(c1.sort.language, c1.sort, c2.sort) is None:
            # i.e. c1 and c2 are disjoint types
            logging.debug('Concept "%s" pruned for type-inconsistency reasons', result)
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

        if isinstance(r, RestrictRole):
            logging.debug('Role "{}" pruned; no direct nesting of restrictions'.format(result))
            return None

        return result

    def create_composition_role(self, r1: Role, r2: Role):

        if r1 == r2:
            return None

        # Compose only on primitives or their inversions
        if (not isinstance(r1, (BasicRole, InverseRole)) or
           not isinstance(r2, (BasicRole, InverseRole))):
            return None

        result = CompositionRole(r1, r2)
        if not self.language.are_vertically_related(r1.sort[1], r2.sort[0]):
            logging.debug('Role "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        return result

    def tests(self, language, interpretations):
        on_r = BasicRole(language.get_predicate("on"))
        ontable_c = BasicConcept(language.get_predicate("ontable"))
        # r = StarRole(BasicRole(language.get_predicate("on")))
        # c = self.create_forall_concept(r, BasicConcept(language.get_predicate("clear")))
        # x = list(interpretations.process_concepts([c]))
        # return x

        # r = StarRole())
        c1 = self.create_exists_concept(on_r, self.top)
        c2 = self.create_not_concept(ontable_c)
        x = list(interpretations.process_concepts([c2, c1]))
        return x


def derive_features(concepts, roles):
    new_features = [BooleanFeature(c) for c in concepts]
    new_features.extend(Numerical1Feature(c) for c in concepts)
    new_features.extend(Numerical2Feature(c1, r, c2) for c1, r, c2 in itertools.product(concepts, roles, concepts))
    return new_features


class InterpretationSet(object):
    def __init__(self, language, states, top, bot):
        self.language = language
        self.states = states
        self.top = top
        self.bot = bot
        self.relevant_predicates = set(p for p in self.language.predicates
                                       if not builtins.is_builtin_predicate(p) and p.arity in (1, 2))
        self.universe = self.compute_universe(states)
        self.cache = self.create_cache_for_samples()

    def generate_extension_trace(self, term):
        assert isinstance(term, (Concept, Role))

        substitution = {}
        trace = bitarray()
        for sid in self.states:
            extension = term.extension(self.cache, sid, substitution)
            self.cache.register_compressed_extension(term, sid, extension)
            trace += extension

        # return create_extension_trace(self.universe, trace, arity)
        return trace

    def build_cache_for_state(self, state, universe):
        cache = dict()
        cache['_atoms_'] = []  # Not sure if we'll need this
        unprocessed = set(self.relevant_predicates)

        for atom in state:
            assert len(atom) in (1, 2, 3)
            name = atom[0]

            if len(atom) == 1:
                cache['_atoms_'].append(name)
            else:
                predicate = self.language.get_predicate(name)
                if predicate in unprocessed:
                    unprocessed.remove(predicate)

                if len(atom) == 2:
                    cache.setdefault(BasicConcept(predicate), set()).add(universe.index(atom[1]))
                else:
                    t = (universe.index(atom[1]), universe.index(atom[2]))
                    cache.setdefault(BasicRole(predicate), set()).add(t)

        cache[self.top] = universe.as_extension()
        cache[self.bot] = set()
        # cache["_index_"] = universe
        cache["_state_"] = state

        # CWA: those predicates not appearing on the state trace are assumed to have empty denotation
        for p in unprocessed:
            term = BasicConcept(p) if p.arity == 1 else BasicRole(p)
            cache[term] = set()

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

    def create_cache_for_samples(self):
        cache = ExtensionCache(self.universe, self.top, self.bot)
        for sid, (_, state) in self.states.items():
            all_terms = self.build_cache_for_state(state, self.universe)
            for term, extension in all_terms.items():
                if isinstance(term, (Concept, Role)):
                    cache.register_extension(term, sid, extension)

        return cache

    def process_term(self, term):
        trace = self.generate_extension_trace(term)
        return self.cache.register_trace(term, trace)

    def process_terms(self, elems):
        return [x for x in elems if self.process_term(x)]


def store_terms(concepts, roles, args):
    os.makedirs('terms', exist_ok=True)
    with open(os.path.join('terms', args.result_filename + '.concepts'), 'w') as f:
        f.write("\n".join(map(str, concepts)))
    with open(os.path.join('terms', args.result_filename + '.roles'), 'w') as f:
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

    # Construct the primitive terms from the input language, and then add the atoms
    concepts, roles = factory.create_primitives()
    concepts = factory.create_atomic_concepts(concepts)
    roles = factory.create_atomic_roles(roles)

    concepts = interpretations.process_terms(concepts)
    roles = interpretations.process_terms(roles)

    # factory.tests(language, interpretations)  # informal tests

    # construct derived concepts and rules obtained with grammar
    c_i, c_j = 0, len(concepts)
    r_i, r_j = 0, len(roles)
    print('\nDeriving concepts and roles using {} iteration(s), starting from {} atomic concepts and {} roles'.
          format(args.k, c_j, r_j))



    for i in range(1, args.k+1):

        old_c, new_c = concepts[0:c_i], concepts[c_i:c_j]
        old_r, new_r = roles[0:r_i], roles[r_i:r_j]

        derived_c = factory.derive_concepts(old_c, new_c, old_r, new_r)
        derived_r = factory.derive_roles(old_c, new_c, old_r, new_r)

        print("it. {}: {} concept(s) and {} role(s) generated".format(i, len(derived_c), len(derived_r)))

        concepts.extend(interpretations.process_terms(derived_c))
        roles.extend(interpretations.process_terms(derived_r))

        c_i, c_j = c_j, len(concepts)
        r_i, r_j = r_j, len(roles)

        print("\t\tof which {} concept(s) and {} role(s) incorporated".format(c_j-c_i, r_j-r_i))

    # profiling.print_snapshot()

    store_terms(concepts, roles, args)
    print('Final output: %d concept(s) and %d role(s)' % (len(concepts), len(roles)))
    print('Number of states in the sample: {}'.format(len(state_samples)))


def configure_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO
    filename = os.path.basename(args.transitions)
    args.result_filename = '.'.join(filename.split('.')[:-1]) + ".{}it".format(args.k)
    filename = os.path.join('logs', args.result_filename + '.log')
    logging.basicConfig(filename=filename, filemode='w', level=level)


def bootstrap(arguments):
    args = parse_arguments(arguments)
    configure_logging(args)
    return args


if __name__ == "__main__":
    _args = bootstrap(sys.argv[1:])
    start = time.process_time()
    # profiling.start()
    main(_args)
    print('CPU time: {:.2f}sec'.format(time.process_time()-start))
    print('Max. memory usage: {:.2f}MB'.format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

