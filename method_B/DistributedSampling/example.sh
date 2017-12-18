#!/bin/bash

SAMPLE_SIZE=7
POPULATION_SIZE=32
BASE_SIZE=4
SEED=1
ITERATIONS=10
PES=2

echo mpirun -n $PES ./optimized/run_experiments --n=$SAMPLE_SIZE --N=$POPULATION_SIZE --k=$BASE_SIZE --seed=$SEED --i=$ITERATIONS
mpirun -n $PES ./optimized/run_experiments --n=$SAMPLE_SIZE --N=$POPULATION_SIZE --k=$BASE_SIZE --seed=$SEED --i=$ITERATIONS
