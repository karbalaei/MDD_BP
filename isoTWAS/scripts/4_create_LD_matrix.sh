#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=run_process_gwas
#SBATCH -c 1
#SBATCH -o logs/4_create_LD_matrix_E.txt
#SBATCH -e logs/4_create_LD_matrix_O.txt
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

Rscript 4_create_LD_matrix.R 

echo "**** Job ends ****"
date