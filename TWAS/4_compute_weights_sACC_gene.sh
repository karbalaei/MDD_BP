#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=twas_compute_weights_sACC_gene_Step4
#SBATCH -c 1
#SBATCH -o logs/4_compute_weights_sACC_gene_O.txt
#SBATCH -e logs/4_compute_weights_sACC_gene_E.txt
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
module load plink/1.90b
module load fusion_twas/github
module load conda_R
#module load gemma


## List current modules
module list

Rscript 4_compute_weights.R -c 10 -r "sACC"



echo "**** Job ends ****"
date