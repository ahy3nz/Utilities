#!/usr/bin/bash
# Given a data file with a list of filenames
# Scp the correct file into the current directory

filename="$1"

while read -r line
do
    scp $RAHMAN:Programs/setup/Bilayer/$line/nvteq*.gro .
done < "$filename"
