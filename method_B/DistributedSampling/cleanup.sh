#!/bin/bash

# Remove libs
scons program=library variant=optimized -j 8

# Remove apps
for program in run_methodD run_methodH run_methodR run_methodP; do 
scons program=$program variant=debug -j 8 -c
if [ "$?" -ne "0" ]; then 
    echo "compile error in $program. exiting."
    exit
fi
scons program=$program variant=optimized -j 8 -c
if [ "$?" -ne "0" ]; then 
    echo "compile error in $program. exiting."
    exit
fi
done

# rm -rf debug
rm -rf optimized

