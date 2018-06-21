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
from tarski.dl import Concept, Role, PrimitiveConcept, PrimitiveRole, InverseRole, StarRole, \
    ConceptCardinalityFeature, MinDistanceFeature, SyntacticFactory, NullaryAtomFeature, NullaryAtom

from extensions import UniverseIndex, ExtensionCache
from parameters import add_domain_goal_parameters
from transitions import read_transitions

signal(SIGPIPE, SIG_DFL)


class TerminologicalFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage, states):
        self.syntax = SyntacticFactory(language)
        self.processor = SemanticProcessor(language, states, self.syntax.top, self.syntax.bot)

    def create_atomic_concepts(self, base_concepts):
        # Start with the base concepts
        concepts = [self.syntax.top, self.syntax.bot] + base_concepts
        concepts = [t for t in concepts if self.processor.process_term(t)]

        # Add Not concepts
        nots = [self.syntax.create_not_concept(c) for c in concepts]
        concepts.extend(t for t in nots if self.processor.process_term(t))
        return concepts

    def create_atomic_roles(self, base_roles):
        # Start with the base roles
        roles = [t for t in base_roles if self.processor.process_term(t)]

        # Add inverse roles
        inverses = [InverseRole(r) for r in roles]
        roles.extend(t for t in inverses if self.processor.process_term(t))

        # Add star roles
        stars = [StarRole(r) for r in roles]
        roles.extend(t for t in stars if self.processor.process_term(t))

        return roles

    def derive_concepts(self, old_c, new_c, old_r, new_r):
        generated = []

        def process(iterator):
            generated.extend(x for x in iterator if self.processor.process_term(x))

        # EXISTS R.C, FORALL R.C
        for roles, concepts in ((new_r, old_c), (old_r, new_c), (new_r, new_c)):
            for fun in (self.syntax.create_exists_concept, self.syntax.create_forall_concept):
                process(fun(r, c) for c, r in itertools.product(concepts, roles))

        # R = R'
        for pairings in [itertools.product(old_r, new_r), itertools.combinations(new_r, 2)]:
            process(self.syntax.create_equal_concept(r1, r2) for r1, r2 in pairings)

        # C AND C'
        for pairings in (itertools.product(new_c, old_c), itertools.combinations(new_c, 2)):
            process(self.syntax.create_and_concept(c1, c2) for c1, c2 in pairings)

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
                process(self.syntax.create_composition_role(r1, r2) for r1, r2 in pairings)

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
        role_restrictions = (self.syntax.create_restrict_role(r, c) for r, c in itertools.product(roles, concepts))
        return (t for t in role_restrictions if self.processor.process_term(t))


class SemanticProcessor(object):
    def __init__(self, language, states, top, bot):
        self.language = language
        self.states = states
        self.top = top
        self.bot = bot
        self.relevant_predicates = set(p for p in self.language.predicates
                                       if not builtins.is_builtin_predicate(p) and p.arity in (0, 1, 2))
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
        unprocessed = set(self.relevant_predicates)

        for atom in state:
            assert len(atom) in (1, 2, 3)
            name = atom[0]

            predicate = self.language.get_predicate(name)
            if predicate in unprocessed:
                unprocessed.remove(predicate)

            if len(atom) == 1:  # i.e. a nullary predicate
                assert predicate.arity == 0
                atom = NullaryAtom(predicate)
                assert atom not in cache
                cache[atom] = True

            elif len(atom) == 2:  # i.e. a unary predicate
                assert predicate.arity == 1
                cache.setdefault(PrimitiveConcept(predicate), set()).add(universe.index(atom[1]))

            else:  # i.e. a binary predicate
                assert predicate.arity == 2
                t = (universe.index(atom[1]), universe.index(atom[2]))
                cache.setdefault(PrimitiveRole(predicate), set()).add(t)

        cache[self.top] = universe.as_extension()
        cache[self.bot] = set()
        # cache["_index_"] = universe
        cache["_state_"] = state

        # CWA: those predicates not appearing on the state trace are assumed to have empty denotation
        for p in unprocessed:
            if p.arity == 1:
                cache[PrimitiveConcept(p)] = set()
            elif p.arity == 2:
                cache[PrimitiveRole(p)] = set()
            else:
                assert p.arity == 0
                cache[NullaryAtom(p)] = False

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
                elif isinstance(term, NullaryAtom):
                    cache.register_nullary_truth_value(term, sid, extension)

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


def create_nullary_features(atoms):
    return [NullaryAtomFeature(a) for a in atoms]


# @profile
def run(config, data):
    assert not data
    problem, language = parse_pddl(config.domain)
    language = add_domain_goal_parameters(problem.domain_name, language)

    logging.info('Loading states and transitions...')
    states, goal_states, transitions = read_transitions(config.sample_file)

    factory = TerminologicalFactory(language, states)

    # Construct the primitive terms from the input language, and then add the atoms
    concepts, roles, atoms = factory.syntax.generate_primitives_from_language()
    logging.info('Primitive (nullary) atoms : {}'.format(", ".join(map(str, atoms))))
    logging.info('Primitive (unary) concepts: {}'.format(", ".join(map(str, concepts))))
    logging.info('Primitive (binary) roles  : {}'.format(", ".join(map(str, roles))))
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
    logging.info('Final output: %d concept(s) and %d role(s)'.format(len(concepts), len(roles)))
    logging.info('Number of states in the sample: {}'.format(len(states)))
    logging.info('Number of concepts with singleton extensions over all states: {}'.format(
        len(factory.processor.singleton_extension_concepts)))

    # Temporarily deactivated, role restrictions very expensive
    # rest = list(factory.create_role_restrictions(concepts, roles))
    # store_role_restrictions(rest, config)
    rest = []
    k = 5
    features = factory.derive_features(concepts, rest, k)
    features += create_nullary_features(atoms)

    logging.info('Final number of features: {}'.format(len(features)))

    # return features, states, goal_states, transitions, factory.processor.cache

    return dict(
        features=features,
        states=states,
        goal_states=goal_states,
        transitions=transitions,
        extensions=factory.processor.cache
    )
