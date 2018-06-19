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
import logging
import itertools
from signal import signal, SIGPIPE, SIG_DFL

import os
import tarski as tsk

from bitarray import bitarray
from tarski.io import FstripsReader
from tarski.syntax import builtins

from extensions import UniverseIndex, ExtensionCache
from parameters import add_domain_goal_parameters
from syntax import Concept, Role, Atom, UniversalConcept, BasicConcept, NotConcept, ExistsConcept, ForallConcept, \
    EqualConcept, BasicRole, InverseRole, StarRole, RestrictRole, ConceptCardinalityFeature, \
    MinDistanceFeature, AndConcept, most_restricted_type, EmptyConcept, CompositionRole, SingletonConcept
from transitions import read_transitions
from util.algorithms import filter_subnodes

signal(SIGPIPE, SIG_DFL)


class TermBox(object):
    """ A simple holder object where we'll classify the different concepts and roles that we generate """

    def __init__(self):
        self.primitive_atoms = []


class TerminologicalFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage):
        self.language = language
        self.universe_sort = language.get_sort('object')
        self.top = UniversalConcept(self.universe_sort.name)
        self.bot = EmptyConcept(self.universe_sort.name)
        self.termbox = TermBox()
        self.processor = None  # Will be set later

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

        for c in self.language.constants():
            concepts.append(SingletonConcept(c.symbol, c.sort))


        atoms = self.termbox.primitive_atoms
        print('{} input (nullary) atoms : {}'.format(len(atoms), list(map(str, atoms))))
        print('{} input (unary) concepts: {}'.format(len(concepts), list(map(str, concepts))))
        print('{} input (binary) roles  : {}'.format(len(roles), list(map(str, roles))))
        return concepts, roles

    # derive new concepts (1-step derivation) from given set of concepts and roles
    def create_atomic_concepts(self, concepts):
        new_concepts = [self.top, self.bot] + concepts + [self.create_not_concept(c) for c in concepts]
        return [term for term in new_concepts if self.processor.process_term(term)]

    def create_atomic_roles(self, roles):
        new_roles = [t for t in roles]

        inverses = [InverseRole(r) for r in roles if not isinstance(r, InverseRole)]

        new_roles.extend(inverses)
        new_roles.extend(StarRole(r) for r in roles if not isinstance(r, StarRole))
        new_roles.extend(StarRole(r) for r in inverses)
        return [term for term in new_roles if self.processor.process_term(term)]

    def derive_concepts(self, old_c, new_c, old_r, new_r):
        generated = []

        def process(iterator):
            generated.extend(x for x in iterator if self.processor.process_term(x))

        # EXISTS R.C, FORALL R.C
        for roles, concepts in ((new_r, old_c), (old_r, new_c), (new_r, new_c)):
            for fun in (self.create_exists_concept, self.create_forall_concept):
                process(fun(r, c) for c, r in itertools.product(concepts, roles))

        # R = R'
        for pairings in [itertools.product(old_r, new_r), itertools.combinations(new_r, 2)]:
            process(self.create_equal_concept(r1, r2) for r1, r2 in pairings)

        # C AND C'
        for pairings in (itertools.product(new_c, old_c), itertools.combinations(new_c, 2)):
            process(self.create_and_concept(c1, c2) for c1, c2 in pairings)

        return generated

    def derive_roles(self, old_c, new_c, old_r, new_r, derive_compositions=True):
        generated = []
        cart_prod = itertools.product
        # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
        # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])

        def process(iterator):
            generated.extend(x for x in iterator if self.processor.process_term(x))

        if derive_compositions:
            for pairings in (cart_prod(old_r, new_r), cart_prod(new_r, new_r)):
                process(self.create_composition_role(r1, r2) for r1, r2 in pairings)

        # for pairings in (cart_prod(new_r, old_c), cart_prod(old_r, new_c), cart_prod(new_r, new_c)):
        #     process(self.create_restrict_role(r, c) for r, c in pairings)

        return generated

    def derive_features(self, concepts, rest, k):
        # new_features = [NonEmptyConceptFeature(c) for c in concepts]
        new_features = []
        new_features.extend(ConceptCardinalityFeature(c) for c in concepts)
        new_features.extend(self.create_distance_features(concepts, rest, k))
        return new_features

    def create_distance_features(self, concepts, rest, k):
        card1_concepts = self.processor.singleton_extension_concepts
        return (MinDistanceFeature(c1, r, c2) for c1, r, c2 in itertools.product(card1_concepts, rest, concepts) if
                c1.depth + r.depth + c2.depth <= k)

    def create_role_restrictions(self, concepts, roles):
        role_restrictions = (self.create_restrict_role(r, c) for r, c in itertools.product(roles, concepts))
        return (t for t in role_restrictions if self.processor.process_term(t))

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

        sort = most_restricted_type(self.language, c1.sort, c2.sort)

        if c1 == c2:
            return None  # No sense in C and C

        if c1 in (self.top, self.bot) or c2 in (self.top, self.bot):
            logging.debug('AND of {} and {} pruned, no sense in AND\'ing with top or bot'.format(c1, c2))
            return None

        if sort is None:
            # i.e. c1 and c2 are disjoint types
            logging.debug('AND of {} and {} pruned for type-inconsistency reasons'.format(c1, c2))
            return None

        return AndConcept(c1, c2, sort)

    def create_equal_concept(self, r1: Role, r2: Role):
        assert isinstance(r1, Role) and isinstance(r2, Role)
        # The sort of the resulting concept will be the most restricted sort between the left sorts of the two roles
        sort = most_restricted_type(self.language, r1.sort[0], r2.sort[0])
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

        # Compose only on primitives or their inversions
        # if (not isinstance(r1, (BasicRole, InverseRole)) or
        #    not isinstance(r2, (BasicRole, InverseRole))):
        #     return None

        result = CompositionRole(r1, r2)

        if not self.language.are_vertically_related(r1.sort[1], r2.sort[0]):
            logging.debug('Role "{}" pruned for type-inconsistency reasons'.format(result))
            return None

        num_comp = len(filter_subnodes(result, CompositionRole))
        if num_comp > 2:
            logging.debug('Role "{}" pruned: number of compositions ({}) exceeds threshold'.format(result, num_comp))
            return None

        return result

    def tests(self, language):
        on_r = BasicRole(language.get_predicate("on"))
        holding_c = BasicConcept(language.get_predicate("holding"))
        ontable_c = BasicConcept(language.get_predicate("ontable"))
        # r = StarRole(BasicRole(language.get_predicate("on")))
        # c = self.create_forall_concept(r, BasicConcept(language.get_predicate("clear")))
        # x = list(interpretations.process_concepts([c]))
        # return x

        # r = StarRole())
        c1 = self.create_exists_concept(on_r, self.top)
        c2 = self.create_not_concept(ontable_c)
        # x = list(interpretations.process_concepts([c2, c1]))

        ##
        not_holding = self.create_not_concept(holding_c)
        a_const = language.get_constant("a")
        a = SingletonConcept("a", a_const.sort)
        not_a = self.create_not_concept(a)

        # Forall(Star(on),Not({a}))
        above = StarRole(on_r)
        below = StarRole(InverseRole(on_r))
        not_above_a = self.create_forall_concept(above, not_a)
        not_below_a = self.create_forall_concept(below, not_a)

        M1 = self.create_and_concept(not_holding, not_a)
        M2 = self.create_and_concept(not_above_a, not_below_a)
        M = self.create_and_concept(M1, M2)



        # (And(Not({a}),Not(holding)), 0)
        # t1 = self.create_and_concept(not_a, not_holding)
        # pt1 = self.processor.process_term(t1)

        z = self.processor.process_term(M)

        return z


class SemanticProcessor(object):
    def __init__(self, language, states, top, bot):
        self.language = language
        self.states = states
        self.top = top
        self.bot = bot
        self.relevant_predicates = set(p for p in self.language.predicates
                                       if not builtins.is_builtin_predicate(p) and p.arity in (1, 2))
        self.universe = self.compute_universe(states)
        self.cache = self.create_cache_for_samples()
        self.singleton_extension_concepts = []

    # @profile
    def generate_extension_trace(self, term):
        assert isinstance(term, (Concept, Role))

        substitution = {}
        trace = bitarray()
        state_extensions = []
        for sid in self.states:
            extension = term.extension(self.cache, sid, substitution)
            trace += extension
            state_extensions.append((sid, extension))

        return trace, state_extensions

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

    @staticmethod
    def compute_universe(states):
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
        if term is None:
            return False
        trace, extensions = self.generate_extension_trace(term)
        if not self.cache.register_trace(term, trace):
            # The trace is equivalent to some other trace already seen, we signal so by returning False
            return False

        singleton_extension_over_all_states = True and term.ARITY == 1
        for sid, ext in extensions:  # We register the compressed individual extensions
            self.cache.register_compressed_extension(term, sid, ext)
            singleton_extension_over_all_states = singleton_extension_over_all_states and ext.count(True) == 1

        if singleton_extension_over_all_states:
            self.singleton_extension_concepts.append(term)
        return True

    def process_terms(self, elems):
        return [x for x in elems if self.process_term(x)]


def store_terms(concepts, roles, args):
    os.makedirs(args.experiment_dir, exist_ok=True)
    with open(os.path.join(args.experiment_dir, 'concepts.txt'), 'w') as f:
        f.write("\n".join(map(str, concepts)))
    with open(os.path.join(args.experiment_dir, 'roles.txt'), 'w') as f:
        f.write("\n".join(map(str, roles)))


def store_role_restrictions(roles, args):
    os.makedirs(args.experiment_dir, exist_ok=True)
    with open(os.path.join(args.experiment_dir, 'role-restrictions.txt'), 'w') as f:
        f.write("\n".join(map(str, roles)))


def parse_pddl(domain_file):
    print('Parsing PDDL domain')
    reader = FstripsReader()
    reader.parse_domain(domain_file)
    problem = reader.problem
    return problem, problem.language


# @profile
def run(config, data):
    assert not data
    problem, language = parse_pddl(config.domain)
    language = add_domain_goal_parameters(problem.domain_name, language)

    print('Loading states and transitions...')
    states, goal_states, transitions = read_transitions(config.sample_file)

    factory = TerminologicalFactory(language)
    factory.processor = SemanticProcessor(language, states, factory.top, factory.bot)

    # Construct the primitive terms from the input language, and then add the atoms
    concepts, roles = factory.create_primitives()
    concepts = factory.create_atomic_concepts(concepts)
    roles = factory.create_atomic_roles(roles)

    # construct derived concepts and rules obtained with grammar
    c_i, c_j = 0, len(concepts)
    r_i, r_j = 0, len(roles)
    print('\nDeriving concepts and roles using {} iteration(s), starting from {} atomic concepts and {} roles'.
          format(config.concept_depth, c_j, r_j))

    for i in range(1, config.concept_depth + 1):
        # Update indexes
        old_c, new_c = concepts[0:c_i], concepts[c_i:c_j]
        old_r, new_r = roles[0:r_i], roles[r_i:r_j]

        print("Starting iteration #{}...".format(i), end='', flush=True)

        derive_compositions = i <= 1
        concepts.extend(factory.derive_concepts(old_c, new_c, old_r, new_r))
        roles.extend(factory.derive_roles(old_c, new_c, old_r, new_r, derive_compositions=derive_compositions))

        c_i, c_j = c_j, len(concepts)
        r_i, r_j = r_j, len(roles)

        print("\t{} new concept(s) and {} new role(s) incorporated".format(c_j-c_i, r_j-r_i))

    # profiling.print_snapshot()

    store_terms(concepts, roles, config)
    print('Final output: %d concept(s) and %d role(s)' % (len(concepts), len(roles)))
    print('Number of states in the sample: {}'.format(len(states)))
    print('Number of concepts with singleton extensions over all states: {}'.format(
        len(factory.processor.singleton_extension_concepts)))
    print('Creating features from {} concepts and {} roles'.format(len(concepts), len(roles)))

    # factory.tests(language)  # some informal tests

    # Temporarily deactivated, role restrictions very expensive
    # rest = list(factory.create_role_restrictions(concepts, roles))
    # store_role_restrictions(rest, config)
    rest = []
    k = 5
    features = factory.derive_features(concepts, rest, k)
    print('Final number of features: {}'.format(len(features)))

    # return features, states, goal_states, transitions, factory.processor.cache

    return dict(
        features=features,
        states=states,
        goal_states=goal_states,
        transitions=transitions,
        extensions=factory.processor.cache
    )
