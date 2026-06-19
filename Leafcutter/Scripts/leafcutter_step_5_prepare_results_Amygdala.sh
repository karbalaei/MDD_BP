#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=DS_leafcutter_Amygdala
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_5_Amygdala.txt
#SBATCH -e ../logs/e_leafcutter_step_5_Amygdala.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


ml leafcutter

./prepare_results.R -o ../leafcutter/covariates_results/step_5_Amygdala.Rdata \
             -m ../data/modSep_Amygdala_for_leafcutter.txt -f 0.05 -c Amygdala  \
             ../leafcutter/covariates/leafcutter_perind_numers.counts.gz \
             ../leafcutter/covariates_results/Amygdala_covariates_cluster_significance.txt \
             ../leafcutter/covariates_results/Amygdala_covariates_effect_sizes.txt \
             ../leafcutter/covariates_results/gencode_hg38
