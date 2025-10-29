#!/bin/bash

## These mkdir steps + ln -s + "mkdir -p logs/Amygdala_gene" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## To avoid having to change the file permissions later
## From https://twitter.com/fellgernon/status/1258455434073124865?s=20
## and https://www.cyberciti.biz/tips/understanding-linux-unix-umask-value-usage.html
umask u=rwx,g=rwx,o= ## equivalent to umask 007

## Required order for running this code:
find ./twas/Amygdala_gene/out_files -type f -name "*.wgt.RDat" > ./twas/Amygdala_gene/Amygdala_list_expression_predict.txt


#Rscript make_score.R [wgt.RDat file] > [SCORE_FILE] 

#plink --bfile [GENOTYPES] --score [SCORE_FILE] 1 2 4


split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.txt ./twas/Amygdala_gene/Amygdala_list_expression_predict.txt ./twas/Amygdala_gene/Amygdala_list_expression_predict_