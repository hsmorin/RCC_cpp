#!/usr/bin/env bash

for i in $(seq $1 $2); do
    cd job$i
    condor_submit job"$i"_row"$3".sub
    cd ..
done
