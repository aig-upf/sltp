import logging
import os

from util.command import execute

_CURRENT_DIR_ = os.path.dirname(os.path.abspath(__file__))


def run(config, data):
    logging.info("Generating SAT encoding problem from {} concept-based features and {} states"
                 .format(len(data.features), len(data.state_ids)))

    command = "{} {} {} {}".format(os.path.join(_CURRENT_DIR_, "encoder", "build_alt_theory"),
                                   config.sat_theory_prefix, config.encoding_k, config.encoding_m)

    execute(command=command.split(' '),
            # stdout=config.sample_file,
            # cwd=config.planner_location
            )
    return dict()
