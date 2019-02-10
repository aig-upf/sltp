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

import os

from .features import parse_all_instances, compute_models
from .util.command import execute
from .returncodes import ExitCode


def run(config, data, rng):
    return extract_features(config, data.sample)


def extract_features(config, sample):
    logging.info("Generating non-redundant concepts from sample set: {}".format(sample.info()))

    parsed_problems = parse_all_instances(config.domain, config.instances)  # Parse all problem instances

    language, nominals, model_cache, infos = compute_models(
        config.domain, sample, parsed_problems, config.parameter_generator)

    all_goal_predicates = set(itertools.chain.from_iterable(info.goal_predicates for info in infos))
    if any(all_goal_predicates != info.goal_predicates for info in infos):
        logging.warning("Not all instances in the training set use the same goal predicate")

    logging.info('Invoking C++ feature generation module'.format())


    # TODO WORK IN PROGRESS

    with open(os.path.join(config.experiment_dir, "problem_metadata.csv")) as f:
        f.write("")

    with open(os.path.join(config.experiment_dir, "primitives.csv")) as f:
        f.write("")

    cmd = os.path.join(config.featuregen_location, "featuregen")
    retcode = execute([cmd, config.experiment_dir])



    types = [s for s in language.sorts if not s.builtin and s != language.Object]
    atoms, concepts, roles = generate_concepts(config, factory, nominals, types, all_goal_predicates)

    logging.info('Final number of features: {}'.format(len(features)))
    # log_concept_denotations(sample.states, concepts, factory.processor.models, config.concept_denotation_filename)

    return ExitCode.Success, dict(
        features=features,
        model_cache=factory.processor.model_cache,
        enforced_feature_idxs=enforced_feature_idxs,
    )



