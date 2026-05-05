#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <system> <n_linear_chains> <config_index>"
    exit 1
fi

system=$1
n_linear_chains=$2
config_index=$3

cd ../${system}_${n_linear_chains}linear/
mkdir -p chain_coordinates/

gcc -O3 -o exe/main_chain_coordinates ../source_code/main_chain_coordinates.c ../source_code/functions.c ../source_code/ran2.c -lm

config_file=$( printf "configurations/config_${system}_${n_linear_chains}linear_%09d.dat" "${config_index}" )
chain_coordinates_file=$( printf "chain_coordinates/chain_coordinates_${system}_${n_linear_chains}linear_%09d.dat" "${config_index}" )

echo "Chain coordinates extraction: configuration ${config_index}" >> log/log_chain_coordinates.txt
{ time ./exe/main_chain_coordinates INPUT_PARAMETERS.DAT $config_file > $chain_coordinates_file ; } >> log/log_chain_coordinates.txt 2>&1 &
