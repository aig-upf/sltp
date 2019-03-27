import logging

from .util.tools import Abstraction
from .features import generate_model_cache
from .returncodes import ExitCode
from .sampling import read_transitions_from_files
from .validator import AbstractionValidator


def run(config, data, rng):
    if config.test_domain is None:
        logging.info("No testing instances were specified")
        return ExitCode.Success, dict()
    logging.info("Testing learnt abstraction on instances: {}".format(config.test_instances))

    abstraction = data.abstraction
    assert isinstance(abstraction, Abstraction)

    sample, _ = read_transitions_from_files(config.test_sample_files)
    model_cache = generate_model_cache(config.test_domain, config.test_instances, sample, config.parameter_generator)

    validator = AbstractionValidator(model_cache, sample, list(sample.expanded))
    flaws = validator.find_flaws(abstraction, 1, check_completeness=False)
    if flaws:
        logging.error("The computed abstraction is not sound & complete".format())
        return ExitCode.AbstractionFailsOnTestInstances, dict()

    logging.info("The computed abstraction is sound & complete in all test instances!".format())
    return ExitCode.Success, dict()
