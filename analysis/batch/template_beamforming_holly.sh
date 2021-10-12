#!/bin/sh

# Name for this job
#$ -N [FILENAME]

#$ -S /bin/sh

# The output from each job - the stuff you would normally see on your screen -
# is written to files. Normally one for error messages, and one for normal output.
# The line below joins the two so everything goes to one file.
#$ -j y

#$ -l vf=[RAM]
#$ -l h_vmem=[RAM]

time /share/apps/matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "cd('[SCRIPTDIR]'); [FUNCTION]('[ARGS]')"
