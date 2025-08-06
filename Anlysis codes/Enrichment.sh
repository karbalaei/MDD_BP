#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=Enrichment
#SBATCH -c 1
#SBATCH -o ../logs/o_Enrichment.txt
#SBATCH -e ../logs/e_Enrichment.txt
#SBATCH --mail-type=ALL

set -e

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
Rscript Enrichment.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.0.0
## available from http://research.libd.org/slurmjobs/
