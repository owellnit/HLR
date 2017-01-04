#!/bin/bash
# Output in "job.out", time and erros in "job.err"
#SBATCH --output=skriptSeq.out
#SBATCH --error=skriptSeq.err
#SBATCH -p west
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
export MPICH_NEMESIS_NETMOD=tcp
fi

mpirun ./partdiff-par 1 1 256 2 1 0.000002332031
