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
import sys

import os
from collections import defaultdict

import tarski as tsk

from bitarray import bitarray
from .models import DLModelFactory, FeatureModel
from sltp.util.misc import compute_universe_from_pddl_model, state_as_atoms
from tarski.syntax.transform.errors import TransformationError
from tarski.syntax.transform.simplifications import transform_to_ground_atoms

from .language import parse_pddl, compute_goal_denotation
from tarski.dl import Concept, Role, InverseRole, StarRole, \
    ConceptCardinalityFeature, MinDistanceFeature, SyntacticFactory, NullaryAtomFeature, compute_dl_vocabulary
from tarski import fstrips

from .extensions import DLDenotationTraceIndex
from .returncodes import ExitCode


class TerminologicalFactory(object):
    def __init__(self, language: tsk.FirstOrderLanguage, processor):
        self.syntax = SyntacticFactory(language)
        self.processor = processor

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

    def derive_concepts(self, old_c, new_c, old_r, new_r, max_size):
        generated = []

        def process(iterator):
            generated.extend(x for x in iterator
                             if x is not None and x.size <= max_size and self.processor.process_term(x))

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

        # NOT C
        # process(self.syntax.create_not_concept(c) for c in new_c)

        return generated

    def derive_roles(self, old_c, new_c, old_r, new_r, max_size, derive_compositions=True):
        generated = []
        cart_prod = itertools.product

        # result.extend([ InverseRole(r) for r in roles if not isinstance(r, InverseRole) ])
        # result.extend([ StarRole(r) for r in roles if not isinstance(r, StarRole) ])

        def process(iterator):
            # -1 because the role will necessarily have to be used within some concept:
            generated.extend(x for x in iterator
                             if x is not None and x.size <= max_size - 1 and self.processor.process_term(x))

        if derive_compositions:
            for pairings in (cart_prod(old_r, new_r), cart_prod(new_r, new_r)):
                process(self.syntax.create_composition_role(r1, r2) for r1, r2 in pairings)

        # for pairings in (cart_prod(new_r, old_c), cart_prod(old_r, new_c), cart_prod(new_r, new_c)):
        #     process(self.create_restrict_role(r, c) for r, c in pairings)

        return generated

    def derive_features(self, concepts, rest, distance_feature_max_complexity):
        # new_features = [NonEmptyConceptFeature(c) for c in concepts]
        feats = []
        feats.extend(ConceptCardinalityFeature(c) for c in concepts)
        k = len(feats)
        logging.info('{} concept cardinality features created'.format(k))

        if distance_feature_max_complexity:
            feats.extend(self.create_distance_features(concepts, rest, distance_feature_max_complexity))
            logging.info('{} min-distance features created'.format(len(feats) - k))

        return feats

    def create_distance_features(self, concepts, rest, k):
        card1_concepts = self.processor.singleton_extension_concepts

        for c1, r, c2 in itertools.product(card1_concepts, rest, concepts):
            if c1.size + r.size + c2.size > k:
                continue
            if c2 in (self.syntax.top, self.syntax.bot):
                continue  # No sense in creating a distance-to-nothing or distance-to-all feature
            if c1 == c2:
                continue

            yield MinDistanceFeature(c1, r, c2)

    def create_role_restrictions(self, concepts, roles):
        all_roles = roles[:]
        role_restrictions = (self.syntax.create_restrict_role(r, c) for r, c in itertools.product(roles, concepts))
        all_roles.extend(t for t in role_restrictions if self.processor.process_term(t))
        return all_roles


class SemanticProcessor(object):
    def __init__(self, states, model_cache):
        self.states = states
        self.model_cache = model_cache
        self.singleton_extension_concepts = []
        self.traces = DLDenotationTraceIndex()

    def process_term(self, term):
        if term is None:
            return False
        trace, all_denotations_are_singleton = self._generate_denotation_trace(term)

        if not self.traces.register_trace(term, trace):
            # The trace is equivalent to some other trace already seen, we signal so by returning False
            return False

        if all_denotations_are_singleton and term.ARITY == 1:
            self.singleton_extension_concepts.append(term)
        return True

    def _generate_denotation_trace(self, term):
        assert isinstance(term, (Concept, Role))
        all_denotations_are_singleton = True
        trace = bitarray()
        for sid in self.states:
            # This doesn't cache the denotation, which is intended, because at this point we still don't know if
            # the term might be redundant, in which case we don't want to clutter the cache with its denotation
            extension = term.denotation(self.model_cache.get_term_model(sid))
            trace += extension
            all_denotations_are_singleton = all_denotations_are_singleton and extension.count(True) == 1
        return trace, all_denotations_are_singleton


def store_terms(concepts, roles, args):
    os.makedirs(args.experiment_dir, exist_ok=True)
    with open(os.path.join(args.experiment_dir, 'concepts.txt'), 'w') as f:
        # f.write("\n".join(map(str, concepts)))
        f.write("\n".join("{} [{}]".format(c, c.size) for c in concepts))
    with open(os.path.join(args.experiment_dir, 'roles.txt'), 'w') as f:
        # f.write("\n".join(map(str, roles)))
        f.write("\n".join("{} [{}]".format(r, r.size) for r in roles))


def store_role_restrictions(roles, args):
    os.makedirs(args.experiment_dir, exist_ok=True)
    with open(os.path.join(args.experiment_dir, 'role-restrictions.txt'), 'w') as f:
        f.write("\n".join(map(str, roles)))


def create_nullary_features(atoms):
    return [NullaryAtomFeature(a) for a in atoms]


# @profile
def collect_all_terms(processor, atoms, concepts, roles):
    """ Collect all simpler terms from the given atoms, concepts and roles """
    assert not atoms  # To be implemented
    all_terms = set()
    for term in itertools.chain(atoms, concepts, roles):
        all_terms.update(term.flatten())

    # Sort elements into categories and from simpler to more complex.
    sorter = lambda elem: elem.size
    all_terms = sorted(all_terms, key=sorter)

    interesting = []
    for t in all_terms:
        try:
            if processor.process_term(t):
                interesting.append(t)
        except KeyError:  # i.e. the denotation of some subterm was not cached
            logging.info('Term "{}" ignored as it involves some empty-denotation subterm'.format(t))
            continue

    def filter_elements(elements, _t):
        return [x for x in elements if isinstance(x, _t)]

    atoms = []  # TODO
    concepts = filter_elements(interesting, Concept)
    roles = filter_elements(interesting, Role)

    return interesting, atoms, concepts, roles


def run(config, data, rng):
    return extract_features(config, data.sample)


def parse_all_instances(domain, instances):
    logging.info("Parsing {} training instances".format(len(instances)))
    return [parse_pddl(domain, instance) for instance in instances]


def compute_instance_information(problem, use_goal_denotation=False):
    static_atoms = defaultdict(list)
    static_predicates = set()
    goal_denotations, goal_predicates = None, set()

    # Compute the universe of each instance - a bit redundant with the above, but should be OK
    universe = compute_universe_from_pddl_model(problem.language)

    # Compute a list with all of the predicate / function symbols from the problem that are static
    # Not sure the fstrips.TaskIndex code is too reliable... static detection seems to be buggy.
    # Let's better do that ourselves with a basic fluent detection routine.
    init_atoms = state_as_atoms(problem.init)
    index = fstrips.TaskIndex(problem.language.name, problem.name)

    index.process_symbols(problem)
    fluent_symbols = {x.head.symbol for x in index.fluent_symbols}

    for atom in init_atoms:
        predicate_name = atom[0]
        if predicate_name not in fluent_symbols:
            static_atoms[predicate_name].append(atom)
            static_predicates.add(predicate_name)

    if use_goal_denotation:
        goal_denotations = defaultdict(list)  # Atoms indexed by their predicate name

        try:
            ground_atoms = transform_to_ground_atoms(problem.goal)
        except TransformationError:
            logging.error("Cannot create goal concepts when problem goal is not a conjunction of ground atoms")
            raise

        for atom in ground_atoms:
            predicate_name = atom[0]
            goal_denotations[predicate_name].append(atom)
            goal_predicates.add(predicate_name)

    return InstanceInformation(universe, static_atoms, static_predicates, goal_denotations, goal_predicates)


class InstanceInformation:
    """ A simple collection of instance data necessary to create the DL models """
    def __init__(self, universe, static_atoms, static_predicates, goal_denotations, goal_predicates):
        self.static_atoms = static_atoms
        self.static_predicates = static_predicates
        self.goal_denotations = goal_denotations
        self.goal_predicates = goal_predicates
        self.universe = universe


def compute_models(domain, sample, parsed_problems, parameter_generator):
    if parameter_generator is not None:
        logging.info('Using user-provided domain parameters and ignoring goal representation')
        use_goal_denotation = False
    else:
        logging.info('Using goal representation and no domain parameters')
        use_goal_denotation = True

    infos = [compute_instance_information(problem, use_goal_denotation) for problem, _, _, _ in parsed_problems]

    # We assume all problems languages are the same and simply pick the first one
    language = parsed_problems[0][0].language
    vocabulary = compute_dl_vocabulary(language)

    nominals, types, model_cache = create_model_cache_from_samples(
        vocabulary, sample, domain, parameter_generator, infos)
    return language, nominals, types, model_cache, infos


def extract_features(config, sample):
    logging.info("Generating concepts and pruning from sample set: {}".format(sample.info()))

    parsed_problems = parse_all_instances(config.domain, config.instances)  # Parse all problem instances

    language, nominals, types, model_cache, infos = compute_models(
        config.domain, sample, parsed_problems, config.parameter_generator)

    processor = SemanticProcessor(sample.states, model_cache)
    factory = TerminologicalFactory(language, processor)

    all_goal_predicates = set(itertools.chain.from_iterable(info.goal_predicates for info in infos))
    if any(all_goal_predicates != info.goal_predicates for info in infos):
        logging.warning("Not all instances in the training set use the same goal predicate")

    if config.concept_generator is not None:
        logging.info('Using set of concepts and roles provided by the user'.format())
        user_atoms, user_concepts, user_roles = config.concept_generator(language)
        _ = collect_all_terms(factory.processor, user_atoms, user_concepts, user_roles)
        # i.e. stick with the user-provided concepts!
        atoms, concepts, roles = user_atoms, user_concepts, user_roles
    elif config.feature_generator is None:
        logging.info('Starting automatic generation of concepts'.format())
        atoms, concepts, roles = generate_concepts(config, factory, nominals, types, all_goal_predicates)
    else:
        atoms, concepts, roles = [], [], []

        # profiling.print_snapshot()

    # store_terms(concepts, roles, config)
    logging.info('Final output: {} concept(s) and {} role(s)'.format(len(concepts), len(roles)))
    logging.info('Number of states in the sample: {}'.format(len(sample.states)))
    logging.info('Number of concepts with singleton extensions over all states: {}'.format(
        len(factory.processor.singleton_extension_concepts)))

    if config.distance_feature_max_complexity:  # If we use Use distance features, we'll want role restrictions
        rest = list(factory.create_role_restrictions(concepts, roles))
        # store_role_restrictions(rest, config)
    else:
        rest = roles

    if config.feature_generator is not None and config.enforce_features:
        raise RuntimeError("Cannot use at the same time options 'feature_generator' and 'enforce_features'")

    enforced_feature_idxs = []
    if config.feature_generator is None:
        features = config.enforce_features(language) if config.enforce_features else []
        enforced_feature_idxs = list(range(0, len(features)))
        features += create_nullary_features(atoms)
        features += factory.derive_features(concepts, rest, config.distance_feature_max_complexity)
    else:
        logging.info('Using user-provided set of features instead of computing them automatically')
        features = config.feature_generator(language)

    logging.info('Final number of features: {}'.format(len(features)))
    # log_concept_denotations(sample.states, concepts, factory.processor.models, config.concept_denotation_filename)

    return ExitCode.Success, dict(
        features=features,
        model_cache=factory.processor.model_cache,
        enforced_feature_idxs=enforced_feature_idxs,
    )


def generate_concepts(config, factory, nominals, types, goal_predicates):
    # Construct the primitive terms from the input language, and then add the atoms
    concepts, roles, atoms = factory.syntax.generate_primitives_from_language(nominals, types, goal_predicates)
    concepts = factory.create_atomic_concepts(concepts)
    roles = factory.create_atomic_roles(roles)
    # construct derived concepts and rules obtained with grammar
    c_i, c_j = 0, len(concepts)
    r_i, r_j = 0, len(roles)
    logging.info('Deriving concepts/roles of max. size {}, starting from {} atomic concepts and {} roles'.
                 format(config.max_concept_size, c_j, r_j))

    i = 1
    max_iterations = config.max_concept_grammar_iterations or sys.maxsize
    while i <= max_iterations:
        # Update indexes
        old_c, new_c = concepts[0:c_i], concepts[c_i:c_j]
        old_r, new_r = roles[0:r_i], roles[r_i:r_j]

        logging.info("Starting iteration #{}...".format(i))

        derive_compositions = i <= 1
        derive_compositions = False  # Temporarily deactivate compositions
        concepts.extend(factory.derive_concepts(old_c, new_c, old_r, new_r, config.max_concept_size))
        roles.extend(factory.derive_roles(old_c, new_c, old_r, new_r, config.max_concept_size,
                                          derive_compositions=derive_compositions))

        c_i, c_j = c_j, len(concepts)
        r_i, r_j = r_j, len(roles)
        num_new_concepts = c_j - c_i
        num_new_roles = r_j - r_i
        logging.info("\t{} new concept(s) and {} new role(s) incorporated".format(num_new_concepts, num_new_roles))
        if num_new_concepts + num_new_roles == 0:
            break
        i += 1
    return atoms, concepts, roles


def log_concept_denotations(states, concepts, models, filename, selected=None):
    selected = selected or concepts
    # selected = ((str(f), f) for f in selected)
    # selected = sorted(selected, key=lambda x: x[0])  # Sort features by name

    with open(filename, 'w') as file:
        for s, concept in itertools.product(states, selected):
            val = models[s].uncompressed_denotation(concept)
            print("s_{}[{}] = {}".format(s, concept, val), file=file)

    logging.info("Concept denotations logged in '{}'".format(filename))


def compute_nominals(domain, parameter_generator):
    # A first parse without all constants to get exactly those constants that we want as nominals
    _, language, nominals, types = parse_pddl(domain)
    if parameter_generator is not None:
        nominals += parameter_generator(language)
    return nominals, types


def create_model_factory(domain, instance, parameter_generator):
    nominals, _ = compute_nominals(domain, parameter_generator)

    # Compute the whole universe of the instance
    problem, language, _, _ = parse_pddl(domain, instance)
    vocabulary = compute_dl_vocabulary(language)
    universe = compute_universe_from_pddl_model(language)

    assert parameter_generator is not None  # Need to refine this to work again when we want goal-denotation DL terms
    goal_denotation = None  # TODO
    # goal_denotation = compute_goal_denotation(problem, use_goal_denotation)
    return problem, DLModelFactory(universe, vocabulary, nominals, goal_denotation)


def create_model_cache_from_samples(vocabulary, sample, domain, parameter_generator, infos):
    """  """
    nominals, types = compute_nominals(domain, parameter_generator)
    model_cache = create_model_cache(vocabulary, sample.states, sample.instance, nominals, infos)
    return nominals, types, model_cache


def create_model_cache(vocabulary, states, state_instances, nominals, infos):
    """ Create a DLModelCache from the given parameters"""
    # First create the model factory corresponding to each instance
    model_factories = []
    for info in infos:
        model_factories.append(DLModelFactory(vocabulary, nominals, info))

    # Then create the model corresponding to each state in the sample
    models = {}
    for sid, state in states.items():
        instance = state_instances[sid]
        models[sid] = model_factories[instance].create_model(state)
    return DLModelCache(models)


class DLModelCache:
    def __init__(self, models):
        """ Create a DLModelCache from a dictionary of precomputed models, indexed by corresponding state """
        self.models = models

    def get_term_model(self, state):
        return self.models[state]

    def get_feature_model(self, state):
        return FeatureModel(self.models[state])
