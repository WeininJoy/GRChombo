#!/bin/bash -l

#SBATCH --nodes 2
## NB cosma7 has 28 cores per node so product of these = 28
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=3
#SBATCH -J ComplexStaticVortex_Weining  #Give it something meaningful.
#SBATCH -o standard_output_file.%J.out
#SBATCH -e standard_error_file.%J.err
#SBATCH -p cosma7-pauper #or some other partition, e.g. cosma, cosma6, etc.
#SBATCH -A dp092 #e.g. dp004
#SBATCH --exclusive
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=r09222023@ntu.edu.tw #PLEASE PUT YOUR EMAIL ADDRESS HERE (without the <>)

module purge
#load the modules used to build your program.
module load intel_comp/2019 
module load intel_mpi/2019 
module load parallel_hdf5/1.10.3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the program
mpirun -np $SLURM_NTASKS ./Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ex params.txt
