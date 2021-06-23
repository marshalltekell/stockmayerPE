#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=cheme          # The account name for the job.
#SBATCH --job-name=mu_0.486	 # The job name.
#SBATCH -c 12                    # The number of cpu cores to use.      
#SBATCH --time=119:59:59          # The time the job will take to run
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL        # send an email when the job starts of ends
#SBATCH --mail-user=mct2180@columbia.edu

module load anaconda/3-2019.03
module load openmpi/gcc/64
mpirun -np 12 --oversubscribe /rigel/home/mct2180/lammps-3Mar20/src/lmp_mpi -in in_RUN.ELEC 
python main_POST.py
# End of script
