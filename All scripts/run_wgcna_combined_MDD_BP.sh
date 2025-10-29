#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=run_wgcna_combined
#SBATCH -c 1
#SBATCH -o ../logs/o_run_wgcna_combined_MDD_BP.txt
#SBATCH -e ../logs/e_run_wgcna_combined_MDD_BP.txt
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command EDIT EDIT
Rscript run_wgcna_combined_MDD_BP.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
