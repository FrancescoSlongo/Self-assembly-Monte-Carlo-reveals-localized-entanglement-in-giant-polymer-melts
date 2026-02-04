#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <system> <n_linear_chains>"
    echo "example: $0 50x50x50 4" 
    exit 1
fi

system=$1
n_linear_chains=$2

INPUT_LATTICE_FILE="../lattice_files/lattice_${system}.dat"
echo "Checking existence of input file $INPUT_LATTICE_FILE"
if [ ! -f "$INPUT_LATTICE_FILE" ] ; then
echo "Fatal error. Input file $INPUT_LATTICE_FILE does not exist."
exit 1
fi


OUTDIR=../${system}_${n_linear_chains}linear
echo "Checking existence of target directory $OUTDIR with input parameters"
if [ ! -d "$OUTDIR" ] ; then
echo "Fatal error. Target directory $OUTDIR does not exist."
exit 1
fi

INPUT_PARAMS_FILE="${OUTDIR}/INPUT_PARAMETERS.DAT"
echo "Checking existence of inpur parameter file $INPUT_PARAMS_FILE"
if [ ! -f "$INPUT_PARAMS_FILE" ] ; then
echo "Fatal error. Input file $INPUT_PARAMS_FILE does not exist."
exit 1
fi


# For good measure, let us (re)create the log and obj subdirectories
mkdir -p ${OUTDIR}/log
mkdir -p ${OUTDIR}/obj
mkdir -p ${OUTDIR}/configurations

cd $OUTDIR

echo "Compiling SAMC executable"
gcc -O3 -o obj/main_sampling ../SAMC_CODE/main.c ../SAMC_CODE/functions.c -lm

echo "Launching SAMC run"
{ time ./obj/main_sampling $INPUT_PARAMS_FILE $INPUT_LATTICE_FILE ; } > log/log_sampling.txt 2>&1
echo "Done"
echo "Output is in directory $OUTDIR/configurations"
