#!/bin/bash -l

#SBATCH --ntasks 4 # The number of cores you need...
#SBATCH -J RealScalarMinkowski_Weining  #Give it something meaningful.
#SBATCH -o standard_output_file.%J.out
#SBATCH -e standard_error_file.%J.err
#SBATCH -p cosma7 #or some other partition, e.g. cosma, cosma6, etc.
#SBATCH -A dp092 #e.g. dp004
#SBATCH --exclusive
#SBATCH -t 2:00:00
#SBATCH --mail-type=END # notifications for job done & fail
#SBATCH --mail-user=r09222023@ntu.edu.tw #PLEASE PUT YOUR EMAIL ADDRESS HERE (without the <>)

module purge
#load the modules used to build your program.
module load intel_comp/2019 
module load intel_mpi/2019 
module load parallel_hdf5/1.10.3


# Run the program
srun Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ex params.txt
