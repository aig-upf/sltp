import logging
from enum import Enum

from .solvers import solve
from .returncodes import ExitCode


class OptimizationPolicy(Enum):
    NUM_FEATURES = 1  # Minimize number of features
    TOTAL_FEATURE_COMPLEXITY = 2  # Minimize the sum of depths of selected features
    NUMERIC_FEATURES_FIRST = 3   # Minimize number of numeric features first, then overall features.


class CompletenessInfo:
    def __init__(self, all_states, all_transitions, optimal_states, optimal_transitions):
        self.all_states = all_states
        self.all_transitions = all_transitions
        self.optimal_states = optimal_states
        self.optimal_transitions = optimal_transitions


def compute_completeness_info(sample, complete_only_wrt_optimal):
    """ Compute optimal states and transitions based on the experiment configuration """
    opt = complete_only_wrt_optimal
    all_states = set(sample.get_sorted_state_ids())
    all_transitions = [(s, s_prime) for s in sample.transitions for s_prime in sample.transitions[s]]

    if opt:
        optimal_states = sample.compute_optimal_states(include_goals=False)
        optimal_transitions = sample.optimal_transitions
    else:
        optimal_states = all_states.copy()
        optimal_transitions = set((x, y) for x, y in all_transitions)  # Consider all transitions as optimal

    cinfo = CompletenessInfo(all_states, all_transitions, optimal_states, optimal_transitions)
    return cinfo


def run_solver(config, data, rng):
    solution = solve(config.experiment_dir, config.cnf_filename, config.maxsat_solver, config.maxsat_timeout)
    if not solution.solved and solution.result == "UNSATISFIABLE":
        return ExitCode.MaxsatModelUnsat, dict()
    else:
        logging.info("MAXSAT solution with cost {} found".format(solution.cost))

    return ExitCode.Success, dict(cnf_solution=solution)
