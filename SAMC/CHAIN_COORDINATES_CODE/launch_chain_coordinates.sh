#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <system> <n_linear_chains> <starting_configuration> <number_of_configurations>"
    echo "example: $0 50x50x50 4 0 10" 
    exit 1
fi

system=$1
n_linear_chains=$2
start_config=$3
n_configs=$4

INPUT_LATTICE_FILE="../lattice_files/lattice_${system}.dat"
echo "Checking existence of input file $INPUT_LATTICE_FILE"
if [ ! -f "$INPUT_LATTICE_FILE" ] ; then
echo "Fatal error. Input file $INPUT_LATTICE_FILE does not exist."
exit 1
fi


OUTDIR=../${system}_${n_linear_chains}linear
echo "Checking existence of target directory $OUTDIR"
if [ ! -d "$OUTDIR" ] ; then
echo "Fatal error. Target directory $OUTDIR does not exist."
exit 1
fi

# For good measure, let us (re)create the chain_coordinates directory
mkdir -p ${OUTDIR}/chain_coordinates

cd $OUTDIR

echo "Compiling CHAIN_COORDINATES executable"
gcc -O3 -o obj/main_chain_coordinates ../CHAIN_COORDINATES_CODE/main.c ../CHAIN_COORDINATES_CODE/functions.c -lm

echo "Launching CHAIN_COORDINATES run"
{ time ./obj/main_chain_coordinates ${system} ${n_linear_chains} ${start_config} ${n_configs} ; } > log/log_chain_coordinates.txt 2>&1 
echo "Done."
echo "Output in directory ${OUTDIR}/chain_coordinates"
