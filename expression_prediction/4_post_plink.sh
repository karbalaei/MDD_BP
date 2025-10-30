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

find ./expression_prediction -type f -name "*.profile" -exec basename {} \; > ./Amygdala_gene_list_prediced.txt.profile

# Define the output directory
OUTPUT_DIR="./expression_prediction"
mkdir -p ${OUTPUT_DIR}

#Rscript make_expression_prediction.R [wgt.RDat file] > [expression_prediction_FILE] 

#plink --bfile [GENOTYPES] --expression_prediction [expression_prediction_FILE] 1 2 4


split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.profile ./Amygdala_gene_list_expression_prediction.profile ./Amygdala_gene_list_expression_prediction_

cd ../sACC_gene

find ./expression_prediction -type f -name "*.wgt.RDat" -exec basename {} \; > ./sACC_gene_list_prediced.txt


#Rscript make_expression_prediction.R [wgt.RDat file] > [expression_prediction_FILE] 

#plink --bfile [GENOTYPES] --expression_prediction [expression_prediction_FILE] 1 2 4


split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.profile ./sACC_gene_list_expression_prediction.profile ./sACC_gene_list_expression_prediction_

# Define the output directory
OUTPUT_DIR="./expression_prediction"
mkdir -p ${OUTPUT_DIR}
