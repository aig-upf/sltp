import logging
import os

from util.command import execute

_CURRENT_DIR_ = os.path.dirname(os.path.abspath(__file__))


def run(config, data):
    state_ids = data.state_ids
    features = data.features
    goal_states = data.goal_states

    logging.info("Generating SAT encoding problem from {} concept-based features and {} states"
                 .format(len(data.features), len(data.state_ids)))

    # string matrix_filename = filename(options.prefix_, K, N, "_matrix.dat", false);
    # string transitions_filename = filename(options.prefix_, K, N, "_transitions.dat", false);
    prefix = os.path.join(config.experiment_dir, "sat")

    command = "{} {} {} {}".format(os.path.join(_CURRENT_DIR_, "encoder", "build_alt_theory"),
                                   prefix, config.encoding_k, config.encoding_m)

    execute(command=command.split(' '),
            # stdout=config.sample_file,
            # cwd=config.planner_location
            )
    return dict()
