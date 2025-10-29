#1_expression_prediction_Amygdala

#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=3_compute_weights_indv_Amygdala_full_gene_1
#SBATCH -c 1
#SBATCH -o logs/predict_genes/compute_weights_indv_Amygdala_full_gene_1_%a_O.txt
#SBATCH -e logs/Amygdala_gene/compute_weights_indv_Amygdala_full_gene_1_%a_E.txt
#SBATCH --mail-type=ALL
#SBATCH --array=1-5000%500

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

echo "SLURM_ARRAY_TASK_ID is ${SLURM_ARRAY_TASK_ID}"


## Load dependencies
module load plink/1.90b
module load fusion_twas/github
module load conda_R


## List current modules
module list


# relative path for FILELIST
FILELIST=$(echo "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/input_ids_1.txt")


## File id and feature name
FEATURENUM=$(awk 'BEGIN {FS="\t"} {print $1}' ${FILELIST} | awk "NR==${SLURM_ARRAY_TASK_ID}")
FEATUREID=$(awk 'BEGIN {FS="\t"} {print $2}' ${FILELIST} | awk "NR==${SLURM_ARRAY_TASK_ID}")

## Define files
FILTBIM="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/bim_files/Amygdala_gene_${FEATURENUM}/filtered_snps_Amygdala_gene_${FEATURENUM}"
TMPFILES="tmp_files/gene_${FEATURENUM}"
OUTFILES="out_files/gene_${FEATURENUM}"
