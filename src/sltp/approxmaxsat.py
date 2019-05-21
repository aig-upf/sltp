import logging

from sltp.errors import CriticalPipelineError
from sltp.util.command import execute

from .returncodes import ExitCode


def run(config, data, rng):
    logging.info("Solving implicit WSAT problem using ad-hoc solver".format())

    # We assume that blai's translator and solver are on the path and called 'wsat-translator' and  'wsat-solver'
    wsat_solver_command = "wsat-solver"
    wsat_translator_command = "wsat-translator"

    # Fill in arguments
    wsat_solver_command += " --implicit"
    if config.complete_only_wrt_optimal: wsat_solver_command += " --complete-only-for-marked-transitions"
    if config.wsat_solver_verbose: wsat_solver_command += " --verbose"

    # Finalize command
    wsat_solver_command += " {}".format(config.sat_theory_prefix)

    # Execute
    retcode = execute(command=wsat_solver_command.split(' '),
                      stdout=config.maxsat_solver_out,
                      )
    if retcode:
        raise CriticalPipelineError("Error running MAXSAT solver. Check the logs. Return code: {}".format(retcode))

    # Parse the output
    with open(config.maxsat_solver_out, 'r') as f:
        out = f.read()

    if "solution: UNSAT" in out:
        return ExitCode.MaxsatModelUnsat, dict()

    # We parse an output like "solution=[0,16,17,2,21]"
    try:
        start = out.find('solution=[') + 10
        end = out.find(']', start)
        features_str = out[start:end]
        features = list(map(int, features_str.split(',')))
        # features_str = re.search('solution=[(.+?)]', out).group(1)
        return ExitCode.Success, dict(selected_feature_idxs=features)
    except (AttributeError, ValueError):
        raise CriticalPipelineError("No solution found in output of MAXSAT solver, in file '{}'".format(
            config.maxsat_solver_out))
