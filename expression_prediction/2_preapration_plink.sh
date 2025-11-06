#!/bin/bash
## These mkdir steps + ln -s + "mkdir -p logs/Amygdala_gene" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## To avoid having to change the file permissions later
## From https://twitter.com/fellgernon/status/1258455434073124865?s=20
## and https://www.cyberciti.biz/tips/understanding-linux-unix-umask-value-usage.html
umask u=rwx,g=rwx,o= ## equivalent to umask 007

## Required order for running this code:
cd ../twas/Amygdala_gene

find ./score -type f -name "*.txt" -exec basename {} \; > ./Amygdala_gene_list_score.txt

# Define the output directory
OUTPUT_DIR="./score"
mkdir -p ${OUTPUT_DIR}

split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.txt ./Amygdala_gene_list_score.txt ./Amygdala_gene_list_score_

cd ../sACC_gene

find ./score -type f -name "*.txt" -exec basename {} \; > ./sACC_gene_list_score.txt

# Define the output directory
OUTPUT_DIR="./score"
mkdir -p ${OUTPUT_DIR}

split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.txt ./sACC_gene_list_score.txt ./sACC_gene_list_score_


