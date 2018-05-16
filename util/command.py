
import os
import subprocess


def run_command(solver, rundir, input_filename):
    error = False
    with open(os.path.join(rundir, 'driver.log'), 'w') as driver_log:
        with open(os.path.join(rundir, 'driver.err'), 'w') as driver_err:
            cmd = solver.command(input_filename)
            print('Executing "{}" on directory "{}"'.format(cmd, rundir))
            retcode = subprocess.call(cmd, cwd=rundir, stdout=driver_log, stderr=driver_err)
    if os.path.getsize(driver_err.name) != 0:
        error = True
    if os.path.getsize(driver_err.name) == 0:  # Delete error log if empty
        os.remove(driver_err.name)
    return error, driver_log.name
