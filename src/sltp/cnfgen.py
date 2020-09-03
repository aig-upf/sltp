import logging
import os
import time

from sltp.separation import extract_features_from_sat_solution, compute_good_transitions, \
    generate_policy_from_sat_solution

from . import SLTP_SRC_DIR
from .util.command import execute, read_file
from .returncodes import ExitCode


def run(config, data, rng):
    from .solvers import solve

    if config.maxsat_encoding != "separation":
        return generate_cnf(config, data)

    good_transitions, features = [], []

    maxiterations = 9999 if config.use_incremental_refinement else 1
    approach_name = "incremental refinement" if config.use_incremental_refinement else "standard"

    solution = None
    for it in range(1, maxiterations+1):
        logging.info(f"Starting iteration {it} of the {approach_name} approach")

        print_good_transition_list(good_transitions, config.good_transitions_filename)
        print_good_features_list(features, config.good_features_filename)

        exitcode, results = generate_cnf(config, data)
        if exitcode == ExitCode.IterativeMaxsatApproachSuccessful:
            # Here we subtract 1 because last iteration was only the check
            print(f"Iterative approach finished successfully after {it-1} iterations")

            # Recover the policy one more time, but doing policy minimization
            _, _, policy = generate_policy_from_sat_solution(config, solution, data.model_cache, minimize_policy=True)

            return ExitCode.Success, dict(transition_classification_policy=policy)

        solution = solve(config.experiment_dir, config.cnf_filename, config.maxsat_solver, config.maxsat_timeout)

        if solution.result == "UNSATISFIABLE":
            return ExitCode.MaxsatModelUnsat, dict()
        elif solution.result == "SATISFIABLE":
            logging.info(f"Possibly suboptimal MAXSAT solution with cost {solution.cost} found")
        else:
            logging.info(f"Optimal MAXSAT solution with cost {solution.cost} found")

        # print_maxsat_solution(solution.assignment, config.wsat_allvars_filename)
        features, good_transitions, policy = generate_policy_from_sat_solution(
            config, solution, data.model_cache, minimize_policy=False)

    return ExitCode.Success, dict(transition_classification_policy=None)


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
    args += ["--use-incremental-refinement"] if config.use_incremental_refinement else []
    args += ["--distinguish-goals"] if config.distinguish_goals else []
    args += ["--cross_instance_constraints"] if config.cross_instance_constraints else []
    retcode = execute([cmd] + args)

    if retcode == 0:
        prepare_maxsat_solver_input(config)

    exitcode = {  # Let's map the numeric code returned by the c++ app into an ExitCode object
        0: ExitCode.Success,
        2: ExitCode.IterativeMaxsatApproachSuccessful
    }.get(retcode, ExitCode.CNFGenerationUnknownError)

    return exitcode, dict()


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
