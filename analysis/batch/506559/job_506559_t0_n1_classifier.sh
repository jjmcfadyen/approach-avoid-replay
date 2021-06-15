#!/bin/bash -l
# Batch script to run a multi-threaded MATLAB job under SGE.
# Request X minutes of wallclock time (format hours:minutes:seconds).  --> FILLED IN BY MATLAB SCRIPT
#$ -l h_rt=0:30:00
# Request X gigabytes of RAM per core. --> FILLED IN BY MATLAB SCRIPT
#$ -l mem=5G
# Request 10 gigabyte of TMPDIR space
#$ -l tmpfs=10G
# Request a number of threads (which will use that number of cores).
# On Myriad you can set the number of threads to a maximum of 36.
#$ -pe smp 12
# Request one MATLAB licence - makes sure your job doesn't start
# running until sufficient licenses are free.
#$ -l matlab=1
# Set the name of the job.  --> FILLED IN BY MATLAB SCRIPT
#$ -N s506559_t0_n1_classifier
# Set the working directory to somewhere in your scratch space.
# This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
# This directory must already exist.  --> FILLED IN BY MATLAB SCRIPT
#$ -wd /home/skgtjm6/Scratch
# Set up directories & modules
module load xorg-utils/X11R7.7
module load matlab/full/r2018b/9.5
module load hdf/5-1.8.15/gnu-4.9.2
module load python2/recommended
module load java/1.8.0_45
module load qt/4.8.6/gnu-4.9.2
module load vtk/6.2.0/gnu-4.9.2
module load fsl/6.0.0
# These echoes output what you are about to run  --> FILLED IN BY MATLAB SCRIPT
/usr/bin/time --verbose matlab -nosplash -nodesktop -nodisplay -r "cd('~/Scratch/2020_RiskyReplay/scripts'); build_classifier('/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/scripts/506559/506559_t0_n1.mat')"
