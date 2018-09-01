import logging
import os

from errors import CriticalPipelineError
from util.command import execute

_CURRENT_DIR_ = os.path.dirname(os.path.abspath(__file__))


def encode(config, data, rng):
    logging.info("Generating SAT encoding problem from {} concept-based features and {} states"
                 .format(len(data.features), len(data.state_ids)))

    command = "{} {} {} {}".format(os.path.join(_CURRENT_DIR_, "encoder", "build_alt_theory"),
                                   config.sat_theory_prefix, config.encoding_k, config.encoding_m)

    retcode = execute(command=command.split(' '),
                      # stdout=config.sample_file,
                      # cwd=config.planner_location
                      )
    if retcode:
        raise CriticalPipelineError("Error running SAT encoder. Check the logs. Return code: {}".format(retcode))
    return dict()


def decode(config, data, rng):
    logging.info("Decoding SAT solution with k={}, m={}".format(config.encoding_k, config.encoding_m))

    command = "{} {} {} {} {}".format(os.path.join(_CURRENT_DIR_, "encoder", "build_alt_theory"), "--decode",
                                      config.sat_theory_prefix, config.encoding_k, config.encoding_m)

    retcode = execute(command=command.split(' '),
                      # stdout=config.sample_file,
                      # cwd=config.planner_location
                      )
    if retcode:
        raise CriticalPipelineError("Error running SAT decoder. Check the logs. Return code: {}".format(retcode))
    return dict()
