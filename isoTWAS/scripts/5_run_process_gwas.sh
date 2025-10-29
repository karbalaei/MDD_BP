#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=run_process_gwas
#SBATCH -c 1
#SBATCH -o logs/5_run_process_gwas_isotwas_E.txt
#SBATCH -e logs/5_run_process_gwas_isotwas_O.txt
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


## Load dependencies
module load conda_R


## List current modules
module list

Rscript 5_process_gwas_isotwas.R 

echo "**** Job ends ****"
date