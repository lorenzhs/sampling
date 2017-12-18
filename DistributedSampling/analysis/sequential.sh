#!/bin/bash

# Init output files
OUTPUT="sequential_comparison"

# Init parameters
PE=1
SAMPLE=$(seq 7 27)
POPULATION=30
BASE=10
ITERATIONS=30

python3 sequential.py -f "output/$OUTPUT" -P $PE -n ${SAMPLE[*]} -N $POPULATION -k $BASE -i $ITERATIONS 
