#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=DS_leafcutter_sACC
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_3_DS_sACC.txt
#SBATCH -e ../logs/e_leafcutter_step_3_DS_sACC.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"



ml leafcutter

leafcutter_ds.R --num_threads 4  ../leafcutter/covariates/leafcutter_perind_numers.counts.gz ../data/modSep_sACC_for_leafcutter.txt -o ../leafcutter/covariates_results/sACC_covariates 

