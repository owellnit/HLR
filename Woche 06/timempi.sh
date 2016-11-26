#!/bin/bash

#SBATCH -o out.txt
#SBATCH -p west
#SBATCH -N 3
#SBATCH --ntasks-per-node=5

/etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
	export MPICH_NEMESIS_NETMOD=tcp
fi

mpiexec ./timempi
