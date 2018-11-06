import os
import subprocess
import sys

from ..solvers import library


class Solution(object):
    def __init__(self):
        self.cost = sys.maxsize
        self.assignment = dict()
        self.solved = False
        self.result = "UNPARSED"

    def parse_result(self, line):
        self.result = line
        self.solved = {"UNSATISFIABLE": False, "UNKNOWN": False, "OPTIMUM FOUND": True}[line]


def parse_maxsat_output(filename):
    solution = Solution()
    with open(filename, 'r') as f:
        # We keep only lines which are the variable assignment (v), the result code (s), or cost bounds (o),
        # the last of which should be the actual cost of the solution, if a solution was found.
        relevant = ((line[0], line[1:].strip()) for line in f if line and line[0] in ('v', 's', 'o'))
        for t, line in relevant:
            if t == 'o':
                solution.cost = min(solution.cost, int(line))
            elif t == 's':
                solution.parse_result(line)
            else:
                assert t == 'v'
                lits = (int(x) for x in line.split(' '))
                solution.assignment = {abs(x): x > 0 for x in lits}

        return solution


def choose_solver(solver):
    if solver == 'wpm3':
        return library.WPM3
    elif solver == 'maxino':
        return library.Maxino
    elif solver == 'openwbo':
        return library.Openwbo
    raise RuntimeError('Unknown solver "{}"'.format(solver))


def solve(rundir, cnf_filename, solver='wpm3'):
    solver = choose_solver(solver)(rundir)
    error, output = run(solver, rundir, cnf_filename, 'maxsat_solver')

    if error:
        raise RuntimeError("There was an error running the MAXSAT solver. Check error logs")

    return parse_maxsat_output(output)


def run(solver, rundir, input_filename, tag=None, stdout=None):
    assert tag or stdout
    error = False
    if stdout is None:
        stdout = os.path.join(rundir, '{}_run.log'.format(tag))

    with open(stdout, 'w') as driver_log:
        with open(os.path.join(rundir, '{}_run.err'.format(tag)), 'w') as driver_err:
            cmd = solver.command(input_filename)
            print('Executing "{}" on directory "{}"'.format(cmd, rundir))
            retcode = subprocess.call(cmd, cwd=rundir, stdout=driver_log, stderr=driver_err)
    if os.path.getsize(driver_err.name) != 0:
        error = True
    if os.path.getsize(driver_err.name) == 0:  # Delete error log if empty
        os.remove(driver_err.name)
    return error, driver_log.name
