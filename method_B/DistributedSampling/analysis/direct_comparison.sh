#!/bin/bash

# Init output files
OUTPUT="direct_comparison"

# Init parameters
PE=4
SAMPLE=(4 5 6 7 8)
POPULATION=32
BASE=10
SEEDS=3

# GNM
python3 direct_comparison.py -f "output/$OUTPUT" -P $PE -n ${SAMPLE[*]} -N $POPULATION -k $BASE -s $SEEDS 
