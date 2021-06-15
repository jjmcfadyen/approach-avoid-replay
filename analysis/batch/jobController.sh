#!/bin/bash
# get list of job files
filelist=job_*.sh
# submit each file using qsub
for file in $filelist
do
   qsub $file
done

