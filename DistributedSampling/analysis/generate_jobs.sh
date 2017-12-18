#!/bin/bash
JOBS="jobs"
RUNS="runs"
LOGS="logs"

# Init parameters
# PE=(1 4 8 16 32 64 128)
PE=(1 4 8)
# SAMPLE=(16 18 20)
SAMPLE=(16)
POPULATION=30
BASE=10
ITERATIONS=30

python generate_jobs.py -P ${PE[*]} -n ${SAMPLE[*]} -N $POPULATION -k $BASE -i $ITERATIONS -jf $JOBS -rf $RUNS -lf $LOGS
