#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=analyze_results_MDD_BP_Michael_method
#SBATCH -c 1
#SBATCH -o logs/o_analyze_results_MDD_BP_Michael_method.txt
#SBATCH -e logs/e_analyze_results_MDD_BP_Michael_method.txt
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
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command EDIT EDIT
Rscript analyze_results_MDD_BP_Michael_method.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
