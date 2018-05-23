
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


def count_file_lines(filename):  # Might be a bit faster with a call to "wc -l"
    i = 0
    with open(filename) as f:
        for i, _ in enumerate(f, 1):
            pass
    return i


def remove_duplicate_lines(filename):
    """ Removes in-place any duplicate line in the file. Will also reorder the lines as a side-effect """
    subprocess.call(['sort', '-u', '-o', filename, filename])


def read_file(filename):
    """ Read a file, line by line, ignoring end-of-line characters"""
    with open(filename) as f:
        for line in f:
            yield line.rstrip('\n')