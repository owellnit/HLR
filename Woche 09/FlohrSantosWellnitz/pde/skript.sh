#!/bin/bash
# Output in "job.out", time and erros in "job.err"
#SBATCH --output=skript.out
#SBATCH --error=skript.err
#SBATCH -p abu
#SBATCH -N 2
#SBATCH --ntasks-per-node=12

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
export MPICH_NEMESIS_NETMOD=tcp
fi

mpirun ./partdiff-par 1 1 256 1 2 512
mpirun ./partdiff-par 1 1 256 2 2 512
mpirun ./partdiff-par 1 2 256 1 2 512
mpirun ./partdiff-par 1 2 256 2 2 512
mpirun ./partdiff-par 1 2 256 2 2 1e-6
