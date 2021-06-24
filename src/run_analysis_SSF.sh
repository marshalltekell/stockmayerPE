#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=cheme          # The account name for the job.
#SBATCH --job-name=SSF	  # The job name.
#SBATCH -c 1                 # The number of cpu cores to use.
#SBATCH --time=23:59:59          # The time the job will take to run
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL        # send an email when the job starts of ends
#SBATCH --mail-user=mct2180@columbia.edu

module load anaconda/3-2019.10
python -u main_SSF.py
# End of script