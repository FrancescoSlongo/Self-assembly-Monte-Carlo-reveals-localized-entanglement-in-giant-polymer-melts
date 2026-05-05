#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <system> <n_linear_chains>"
    exit 1
fi

system=$1
n_linear_chains=$2

cd ../${system}_${n_linear_chains}linear/

gcc -O3 -o exe/main_sampling ../source_code/main_sampling.c ../source_code/functions.c ../source_code/ran2.c -lm

{ time ./exe/main_sampling INPUT_PARAMETERS.DAT ; } > log/log_sampling.txt 2>&1 &
