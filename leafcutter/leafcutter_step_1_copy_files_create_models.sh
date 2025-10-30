#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=leafcutter_step_1_copy_files_create_models
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_1_copy_files_create_models.txt
#SBATCH -e ../logs/e_leafcutter_step_1_copy_files_create_models.txt
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript leafcutter_step_1_copy_files_create_models.R

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/leafcutter/junc_files

find "$(pwd)" -name "*.bed" > leafcutter_DS_juncfiles_covariates.txt


echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.0.0
## available from http://research.libd.org/slurmjobs/
