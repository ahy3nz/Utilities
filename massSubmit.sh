#!/usr/bin/bash
# Given a data file with a list of filenames
# CD into that directory and submit the appropriate pbs script
# requires the filename as argument 1 and the extnesion (ST or MD) as argument 2

filename="$1"
extension="$2"

first=true
while read -r line
do
    cd ~/Trajectories/$line
    if [ "$first" = true ];
    then
        echo "$line"
        jobid=$(qsub "$line""$extension"pbs.pbs)
        first=false
    else
        echo "$line"
        jobid=$(qsub -W depend=afterany:$jobid "$line""$extension"pbs.pbs)
    fi
done < "$filename"
