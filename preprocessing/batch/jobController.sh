#!/bin/bash
# get list of job files
filelist=~/Scratch/2021_HrvojeProject/batch/job_*.sh
# submit each file using qsub
for file in $filelist
do
   qsub $file
done
