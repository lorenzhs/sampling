#!/bin/bash

# Build libs
scons program=library variant=optimized -j 8

# Build apps
for program in run_methodH run_methodR run_methodSH; do
# for program in run_experiments; do
scons program=$program variant=debug -j 8
if [ "$?" -ne "0" ]; then
     echo "compile error in $program. exiting."
     exit
fi
scons program=$program variant=optimized -j 8
if [ "$?" -ne "0" ]; then
    echo "compile error in $program. exiting."
    exit
fi
done
