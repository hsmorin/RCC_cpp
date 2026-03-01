#!/bin/bash

i=0
numJobs=$1

while [ $i -lt $numJobs ]; do
    python genDF.py ${i} L $2
    ((i++))
done
