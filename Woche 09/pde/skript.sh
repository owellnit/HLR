#!/bin/bash
# Output in "job.out", time and erros in "job.err"
#SBATCH --output=skript.out
#SBATCH --error=skript.err
#SBATCH -p west
#SBATCH -N 2
#SBATCH --ntasks-per-node=12

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
export MPICH_NEMESIS_NETMOD=tcp
fi

mpirun ./partdiff-par 1 1 128 2 2 128
