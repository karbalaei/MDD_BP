#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=Compare_t_BP_t_MDDBP_t_MDD
#SBATCH -c 1
#SBATCH -o ../logs/o_Compare_t_BP_t_MDDBP_t_MDD.txt
#SBATCH -e ../logs/e_Compare_t_BP_t_MDDBP_t_MDD.txt
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
Rscript Compare_t_BP_t_MDDBP_t_MDD.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.0.0
## available from http://research.libd.org/slurmjobs/
