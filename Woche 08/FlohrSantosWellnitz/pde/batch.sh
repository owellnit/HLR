#!/bin/bash
# Use "west" partition.
#SBATCH --partition=west
# Output in "job.out", time and erros in "job.err"
#SBATCH --output=job.out
#SBATCH --error=job.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
export MPICH_NEMESIS_NETMOD=tcp
fi

mpirun -np 4 ./partdiff-par 12 2 512 2 2 1024
