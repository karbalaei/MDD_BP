#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=isotwas_Amygdala
#SBATCH -c 1
#SBATCH -o logs/Amygdala/3_compute_isotwas_Amygdala_%a_O.txt
#SBATCH -e logs/Amygdala/3_compute_isotwas_Amygdala_%a_E.txt
#SBATCH --mail-type=ALL
#SBATCH --array=1-6499%500
#SBATCH --time=3-00:00:00

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

echo "SLURM_ARRAY_TASK_ID is ${SLURM_ARRAY_TASK_ID}"
## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load plink/1.90b
module load conda_R

## List current modules for reproducibility
module list

# relative path for FILELIST
gene_list_isotwas="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/Amygdala_transcript/gene_list_isotwas_Amygdala.txt"

## File id and feature name
gene_index=$(sed 1d ${gene_list_isotwas} | awk 'BEGIN {FS="\t"} {print $1}' ${gene_list_isotwas} | awk "NR==${SLURM_ARRAY_TASK_ID}")
gene_id=$(sed 1d ${gene_list_isotwas} |awk 'BEGIN {FS="\t"} {print $2}' ${gene_list_isotwas} | awk "NR==${SLURM_ARRAY_TASK_ID}")

#debugging command, check what the names are.
echo "gene_index: $gene_index"
echo "gene_id: $gene_id"

## Edit with your job command EDIT EDIT
Rscript Compute_isotwas.R \
	--region Amygdala \
	--gene_index ${gene_index} \
	--gene_id "${gene_id}"

echo "**** Job ends ****"
date

