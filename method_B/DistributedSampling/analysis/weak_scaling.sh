#!/bin/bash

# Init output files
OUTPUT="weak_scaling"

# Init parameters
PE=(1 2 3 4)
SAMPLE=100000000
POPULATION=52
BASE=10
SEEDS=3

# TODO Bad weak scaling comes from base case sampling
python3 weak_scaling.py -f "output/$OUTPUT" -P ${PE[*]} -n $SAMPLE -N $POPULATION -k $BASE -s $SEEDS 
