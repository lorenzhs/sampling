#!/bin/bash

# Init output files
OUTPUT="parallel_comparison"

# Init parameters
PE=(1 2 4 8 16 32 64 128)
SAMPLE=(16 18 20)
POPULATION=30
BASE=10
ITERATIONS=30

python parallel.py -f "output/$OUTPUT" -P ${PE[*]} -n ${SAMPLE[*]} -N $POPULATION -k $BASE -i $ITERATIONS 
