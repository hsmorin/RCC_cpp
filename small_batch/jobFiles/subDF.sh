#!/usr/bin/env bash

for i in $(seq $1 $2); do
    cd job$i
    condor_submit_dag job$i.dag
    cd ..
done
