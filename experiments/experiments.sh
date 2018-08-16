#!/usr/bin/env bash -ve

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
### Set memory limit.
#SBATCH --mem-per-cpu=10G
### Number of tasks.
#SBATCH --array=1-4
### Adjustment to priority ([-2147483645, 2147483645]).
#SBATCH --nice=2000
### Send mail? Mail type can be e.g. NONE, END, FAIL, ARRAY_TASKS.
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=guillem.frances@unibas.ch


module purge
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no module load Python/3.5.2-goolf-1.7.20
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no module load GCC/5.4.0-2.26
source ${HOME}/lib/virtualenvs/concepts/bin/activate

./cluster.py --exp experiments --task ${SLURM_ARRAY_TASK_ID}
