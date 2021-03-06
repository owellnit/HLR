#!/bin/bash

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=messung.out
# Use "west" partition.
#SBATCH --partition=west

run_output="messung.log"
TIMEFORMAT='%E'

for (( i = 1; i <= 12; i++ )); do
	echo "$i Thread(s):" >> $run_output
	echo "$i Thread(s): "

	for (( j = 1; j <= 3; j++ )); do
		time_for_run=$(time (./partdiff-posix $i 2 512 2 2 1024 >> $run_output 2>&1) 2>&1)
		echo Lauf $j: $time_for_run seconds
	done
done
