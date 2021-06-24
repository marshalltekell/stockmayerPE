#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=cheme          # The account name for the job.
#SBATCH --job-name=mu_0.347	 # The job name.
#SBATCH -c 12                    # The number of cpu cores to use.      
#SBATCH --time=119:59:59          # The time the job will take to run
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL        # send an email when the job starts of ends
#SBATCH --mail-user=mct2180@columbia.edu

module load anaconda/3-2019.10
module load openmpi/gcc/64/4.0.0
python main_PRE.py
mpirun -np 12 --oversubscribe /moto/home/mct2180/mylammps/src/lmp_mpi -in in_EQUIL.ONEAN
# End of script
