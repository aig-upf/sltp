#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import yaml


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp', type=str,  required=True, help="Generate the cluster bash script to be run with squeue "
                                                 "from the given experiment config file")
    parser.add_argument('--task', type=int, default=None, help="Task ID, if in run mode")

    return parser


def main(parser, args):
    args = parser.parse_args(args)
    with open("{}.yml".format(args.exp), 'r') as f:
        data = yaml.load(f)
        experiments = data["configurations"]

    if not experiments:
        raise RuntimeError("No experiments found in the configuration file")

    if args.task is None:
        generate_script(time=60, mem=64, num_tasks=len(experiments), experiment_set=args.exp)
    else:
        # Simply run the whole thing!
        if args.task - 1 > len(experiments):
            raise RuntimeError("Task ID #{} not defined on experiment set {}.".format(args.task, args.exp))

        import runner
        d, e = experiments[args.task].split(" ")
        runner.run(["-d", d, "-e", e])


def generate_script(num_tasks, time, mem, experiment_set):
    tpl = """#!/usr/bin/env bash -ve

### Set name.
#SBATCH --job-name=gen-feats
### Redirect stdout and stderr.
#SBATCH --output=slurm.log
#SBATCH --error=slurm.err
### Let later steps append their logs to the output and error files.
#SBATCH --open-mode=append
### Set partition.
#SBATCH --partition=infai
### Set quality-of-service group.
#SBATCH --qos=infai
### Set time limit (in min).
#SBATCH --time={time}
### Set memory limit.
#SBATCH --mem-per-cpu={mem}G
### Number of tasks.
#SBATCH --array=1-{num_tasks}
### Adjustment to priority ([-2147483645, 2147483645]).
#SBATCH --nice=2000
### Send mail? Mail type can be e.g. NONE, END, FAIL, ARRAY_TASKS.
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=guillem.frances@unibas.ch


module purge
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no module load Python/3.5.2-goolf-1.7.20
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no module load GCC/5.4.0-2.26
source ${{HOME}}/lib/virtualenvs/concepts/bin/activate

./cluster.py --exp {experiment_set} --task ${{SLURM_ARRAY_TASK_ID}} > output_${{SLURM_ARRAY_TASK_ID}}
"""
    filename = "{}.sh".format(experiment_set)
    with open(filename, "w") as f:
        f.write(tpl.format(time=time, mem=mem, num_tasks=num_tasks, experiment_set=experiment_set))
    print("Written cluster script {}".format(filename))


if __name__ == "__main__":
    main(create_parser(), sys.argv[1:])
