#!/bin/bash

FILES="jobs/job_*"
for f in $FILES 
do
    echo $f
    # llsubmit $f
done
