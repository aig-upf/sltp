import logging
import os

from errors import CriticalPipelineError
from solvers import library, common

_CURRENT_DIR_ = os.path.dirname(os.path.abspath(__file__))


def choose_solver(solver):
    if solver == 'glucose':
        return library.Glucose
    elif solver == 'glucose-syrup':
        return library.GlucoseSyrup
    elif solver == 'openwbo':
        return library.Openwbo
    raise RuntimeError('Unknown solver "{}"'.format(solver))


def run(config, data, rng):
    # solution = solve(config.experiment_dir, config.cnf_filename, 'wpm3')
    # solution = solve(config.experiment_dir, config.cnf_filename, 'maxino')

    solver = "glucose-syrup"
    # solver = "openwbo"

    solver = choose_solver(solver)
    error, output = common.run(solver, config.experiment_dir, config.sat_theory_filename,
                               stdout=config.sat_solution_filename)

    if error:
        raise RuntimeError("There was an error running the SAT solver. Check error logs")

    solution = common.parse_maxsat_output(output)  # Sat output just happens to be a subset of maxsat output

    if not solution.solved and solution.result == "UNSATISFIABLE":
        raise CriticalPipelineError("SAT encoding is UNSATISFIABLE")
    else:
        logging.info("SAT solution found!".format())

    return dict()
