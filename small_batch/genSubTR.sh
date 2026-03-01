#!/bin/bash

i=0
numJobs=$1

while [ $i -lt $numJobs ]; do
    python genSubTR.py ${i} $2 L
    ((i++))
done
