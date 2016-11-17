#!/bin/bash

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=messung_mutex.out
#SBATCH --partition=west

run_output="messung_mutex.log"
TIMEFORMAT='%E'

for (( j = 1; j <= 3; j++ )); do
    time_for_run=$(time (./04-PDE/partdiff-seq 1 2 512 2 2 1024 >> $run_output 2>&1) 2>&1)
    echo Lauf $j: $time_for_run seconds
done

for (( i = 1; i <= 12; i++ )); do
	echo "$i Thread(s):" >> $run_output
	echo "$i Thread(s): "

	for (( j = 1; j <= 3; j++ )); do
		time_for_run=$(time (./pde_mutex/partdiff-posix $i 2 512 2 2 1024 >> $run_output 2>&1) 2>&1)
		echo Lauf $j: $time_for_run seconds
	done
done

