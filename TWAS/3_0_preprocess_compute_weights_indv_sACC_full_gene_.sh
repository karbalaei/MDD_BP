#!/bin/bash

## These mkdir steps + ln -s + "mkdir -p logs/sACC_gene" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## To avoid having to change the file permissions later
## From https://twitter.com/fellgernon/status/1258455434073124865?s=20
## and https://www.cyberciti.biz/tips/understanding-linux-unix-umask-value-usage.html
umask u=rwx,g=rwx,o= ## equivalent to umask 007

## Required order for running this code:
cd sACC_gene

split -l 4500 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.txt input_ids.txt input_ids_


## Create symbolic link which is necessary to have accompany between TWAS-FUSION r script and gemma


ln -s  /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/sACC_gene  /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/sACC_gene/sACC_gene
ln -s  /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/sACC_gene  /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/sACC_gene/output


