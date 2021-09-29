#!/bin/bash

##############################################################

# Script Name:        INSIDER_Pipeline.sh
# Author:             Aidan Tay

# Description: This is the general workflow of INSIDER

################### Workspace & Notes #########################

## sh INSIDER_Pipeline.sh \
##     ../../../data/FMDV_serotype_O_full_23062021_subset.fasta \
##     3 \
##     test

################### Dependencies ##############################

################### Global Variables ##########################

## Read command line arguments. Arguments SHOULD be the
## input file and the output directory (among other files)
FASTA_FILE=$1
KMER_LEN=$2
OUTPUT_DIR=$3

## Output directories
KMER_FREQ_DIR=freqs
CID_FILE=cId
RESULTS_DIR=results
CURRDATE="$(date +'%Y-%m-%d %H:%M:%S')"

################### Functions #################################

#################### Main #####################################

## Zip up the src files we need
zip -q \
    -x "**/__pycache__/*" \
    -r src/modules.zip src

echo '###############################################################'
echo
echo 'This analysis was run on:'
echo ${CURRDATE}
echo
echo '###############################################################'
echo
echo 'PWD:'
echo $PWD
echo
echo '###############################################################'
echo
echo 'Args:'
echo $FASTA_FILE
echo $KMER_LEN
echo $OUTPUT_DIR
echo
echo '###############################################################'
echo

## Get the oligonucleotide frequencies for each sequence
time python bin/calculate_kmer_frequencies.py \
    split \
    -f ${FASTA_FILE} \
    -k ${KMER_LEN} \
    -n \
    -o ${OUTPUT_DIR}/${KMER_FREQ_DIR}

echo
echo '###############################################################'
echo

## Cluster related genomic signatures
time python insider_cluster.py \
    consensus \
    --freqDir ${OUTPUT_DIR}/${KMER_FREQ_DIR} \
    -o ${OUTPUT_DIR}/${CID_FILE}.txt

echo
echo '###############################################################'
echo

## Analyse genomic signatures
time python insider_analyse.py \
    main \
    --freqDir ${OUTPUT_DIR}/${KMER_FREQ_DIR} \
    --cIdFile ${OUTPUT_DIR}/${CID_FILE}.txt \
    -o ${OUTPUT_DIR}/${RESULTS_DIR}.txt
