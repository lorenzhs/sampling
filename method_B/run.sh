#!/bin/bash

universe=$((2**30))

SUFF=
if [[ -f rand32-pgo && -f rand32S-pgo && -f rand64-pgo && -f rand64S-pgo ]]; then
    echo "Using PGO'd binaries"
    SUFF=-pgo
fi

# sequential
for logsamples in {7..27}; do
    echo "==============================================="
    samples=$((2**$logsamples))
    iterations=$((2**30 / $samples))
    echo "Running with 2^$logsamples samples ($iterations iterations)"
    time ./rand32$SUFF -q -t 1 -k $samples -n $universe -i $iterations $@
    time ./rand32S$SUFF -q -t 1 -k $samples -n $universe -i $iterations $@
    echo "==============================================="
    echo ""
done

echo ""
echo "================================================================"
echo "================================================================"
echo "======================= Now with N=2^50 ========================"
echo "================================================================"
echo "================================================================"
echo ""
echo ""
universe=$((2**50))

# sequential
for logsamples in {7..27}; do
    echo "==============================================="
    samples=$((2**$logsamples))
    iterations=$((2**30 / $samples))
    echo "Running with 2^$logsamples samples ($iterations iterations)"
    time ./rand64$SUFF -q -t 1 -k $samples -n $universe -i $iterations $@
    time ./rand64S$SUFF -q -t 1 -k $samples -n $universe -i $iterations $@
    echo "==============================================="
    echo ""
done
