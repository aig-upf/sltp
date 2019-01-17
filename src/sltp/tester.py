#!/usr/bin/env python3
import logging

from .features import parse_all_instances, compute_models
from .returncodes import ExitCode
from .sampling import read_transitions_from_files
from .validator import AbstractionValidator


def run(config, data, rng):
    abstraction = {"abstract_actions": data.abstract_actions,
                   "selected_features": data.selected_features,
                   "features": data.features}

    if config.test_domain is None:
        logging.info("No testing instances were specified")
        return ExitCode.Success, dict()

    logging.info("Testing learnt abstraction on instances: ".format(config.test_instances))
    sample, goals_by_instance = read_transitions_from_files(config.test_sample_files)

    parsed_problems = parse_all_instances(config.test_domain, config.test_instances)
    language, nominals, types, model_cache, infos = compute_models(
        config.domain, sample, parsed_problems, config.parameter_generator)

    # we don't care about the order of validation
    validator = AbstractionValidator(model_cache, sample, list(sample.expanded))
    flaws = validator.find_flaws(abstraction, 1, check_completeness=False)
    if flaws:
        logging.error("The computed abstraction is not sound & complete".format())
        return ExitCode.AbstractionFailsOnTestInstances, dict()

    logging.info("The computed abstraction is sound & complete in all test instances!".format())
    return ExitCode.Success, dict()
