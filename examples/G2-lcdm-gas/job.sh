#!/bin/bash 
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=harry.potter@hogwarts.edu
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name G2-lcdm-gas

echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

mpiexec -np $SLURM_NPROCS  ./Gadget4 param.txt

