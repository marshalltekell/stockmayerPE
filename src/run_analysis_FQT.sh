#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=cheme          # The account name for the job.
#SBATCH --job-name=fqt	  # The job name.
#SBATCH -c 6	                 # The number of cpu cores to use.
#SBATCH --time=11:59:59          # The time the job will take to run
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL        # send an email when the job starts of ends
#SBATCH --mail-user=mct2180@columbia.edu

module load anaconda/3-2019.10
module load openmpi/gcc/64/4.0.0
export CFLAGS="-I/moto/opt/anaconda3-2019.10/lib/python3.7/site-packages/numpy/core/include $CFLAGS"
swig -python cfunctions.i
python setup.py build_ext --inplace
mpirun -np 6 --oversubscribe python -u main_FQT.py
# End of script