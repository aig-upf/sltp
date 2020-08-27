import logging
import os
import time

from sltp.separation import extract_features_from_sat_solution, compute_good_transitions
from sltp.util.naming import compute_info_filename

from . import SLTP_SRC_DIR
from .util.command import execute, read_file
from .returncodes import ExitCode


def run(config, data, rng):
    from .solvers import solve

    if config.maxsat_encoding != "separation":
        return generate_cnf(config, data)

    good_transitions, good_features = [], []

    it = 1
    while True:
        print_good_transition_list(good_transitions, config.good_transitions_filename)
        print_good_features_list(good_features, config.good_features_filename)

        exitcode, data = generate_cnf(config, data)
        if exitcode == ExitCode.IterativeMaxsatApproachSuccessful:
            print(f"Iterative approach finished successfully after {it} iterations")

        if exitcode != ExitCode.Success:
            return exitcode, data

        solution = solve(config.experiment_dir, config.cnf_filename, config.maxsat_solver, config.maxsat_timeout)

        if solution.result == "UNSATISFIABLE":
            return ExitCode.MaxsatModelUnsat, dict()

        logging.info(f"MAXSAT solution with cost {solution.cost} found")
        # print_maxsat_solution(solution.assignment, config.wsat_allvars_filename)
        good_transitions = compute_good_transitions(solution.assignment, config.wsat_varmap_filename)
        good_features = extract_features_from_sat_solution(config, solution)

        it += 1


def print_good_transition_list(good_txs, filename):
    with open(filename, 'w') as f:
        for s, sprime in good_txs:
            print(f"{s} {sprime}", file=f)


def print_good_features_list(good_features, filename):
    with open(filename, 'w') as f:
        print(' '.join(str(f.id) for f in good_features), file=f)


def generate_cnf(config, data):
    # Invoke C++ feature generation module
    logging.info('Invoking C++ CNF generation module'.format())
    featuregen_location = os.path.join(SLTP_SRC_DIR, "..", "features")
    cmd = os.path.realpath(os.path.join(featuregen_location, "cnfgen"))
    args = ["--workspace", config.experiment_dir]
    args += ["--enforce-features", ",".join(map(str, data.in_goal_features))] if data.in_goal_features else []
    args += ["--prune-redundant-states"] if config.prune_redundant_states else []
    args += ["--encoding", config.maxsat_encoding]
    args += ["--distinguish-transitions-locally"] if config.distinguish_transitions_locally else []
    args += ["--use-equivalence-classes"] if config.use_equivalence_classes else []
    args += ["--use-feature-dominance"] if config.use_feature_dominance else []
    args += ["--v_slack", str(config.v_slack)]
    retcode = execute([cmd] + args)
    if retcode != 0:
        return ExitCode.CNFGenerationUnknownError, dict()
    prepare_maxsat_solver_input(config)
    return ExitCode.Success, dict()


def prepare_maxsat_solver_input(config):
    """ Read off the output of the C++ CNF generation module and do some transformations, mainly replacing the string
    "TOP" with the actual maxsat TOP integer value.
    """
    with open(config.top_filename, "r") as f:
        top, nvars, nclauses = [int(x) for x in f.read().strip("\n").split(' ')]
    with open(config.cnf_filename, "w") as output:
        print("c WCNF model generated on {}".format(time.strftime("%Y%m%d %H:%M:%S", time.localtime())), file=output)
        print("c Next line encodes: wcnf <nvars> <nclauses> <top>", file=output)
        # p wcnf nvars nclauses top
        print(f"p wcnf {nvars} {nclauses} {top}", file=output)
        for line in read_file(config.cnf_filename + ".tmp"):
            print(line.replace("TOP", str(top)), file=output)
