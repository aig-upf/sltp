import logging
import os
import time

from . import SLTP_SRC_DIR
from .util.command import execute, read_file
from .returncodes import ExitCode


def run(config, data, rng):

    # Invoke C++ feature generation module
    logging.info('Invoking C++ CNF generation module'.format())
    featuregen_location = os.path.join(SLTP_SRC_DIR, "..", "features")
    cmd = os.path.realpath(os.path.join(featuregen_location, "cnfgen"))
    args = ["--workspace", config.experiment_dir]

    if data.in_goal_features:
        args += ["--enforce-features", ",".join(map(str, data.in_goal_features))]

    if config.prune_redundant_states:
        args.append("--prune-redundant-states")

    args += ["--encoding", config.maxsat_encoding]

    args += ["--distinguish-transitions-locally"] if config.distinguish_transitions_locally else []
    args += ["--use-equivalence-classes"] if config.use_equivalence_classes else []
    args += ["--use-feature-dominance"] if config.use_feature_dominance else []
    args += ["--v_slack", str(config.v_slack)]

    retcode = execute([cmd] + args)

    if retcode != 0:
        return ExitCode.CNFGenerationUnknownError, dict()

    # Read off the output of the module and do some transformations
    tmpfilename = config.cnf_filename + ".tmp"

    with open(config.top_filename, "r") as f:
        top, nvars, nclauses = [int(x) for x in f.read().strip("\n").split(' ')]

    with open(config.cnf_filename, "w") as output:
        print("c WCNF model generated on {}".format(time.strftime("%Y%m%d %H:%M:%S", time.localtime())), file=output)
        print("c Next line encodes: wcnf <nvars> <nclauses> <top>".format(nvars, nclauses, top), file=output)
        # p wcnf nvars nclauses top
        print("p wcnf {} {} {}".format(nvars, nclauses, top), file=output)
        for line in read_file(tmpfilename):
            print(line.replace("TOP", str(top)), file=output)

    return ExitCode.Success, dict()
