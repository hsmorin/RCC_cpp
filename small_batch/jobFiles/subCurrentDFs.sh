#!/usr/bin/env bash

for jobInd in "$@"; do
    cd job"$jobInd"
    condor_submit_dag job"$jobInd".dag
    cd ..
done
