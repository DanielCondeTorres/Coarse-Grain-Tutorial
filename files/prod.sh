#!/bin/bash
#SBATCH -t  00:30:00 # execution time. Ex: 1/2 hour
#SBATCH --mem-per-cpu=1GB
#SBATCH -n 1 -c 1# number of tasks, number of cores
#SBATCH --ntasks-per-node=1
module load gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
srun gmx_mpi mdrun -pin on -cpi -noappend -s prod.tpr 
