#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=1G
#SBATCH --job-name=compute_predict_expression_Amygdala_full_gene_1
#SBATCH -c 1
#SBATCH -o ../twas/logs/predict_expression/plink/compute_predict_expression_Amygdala_full_gene_1_%a_O.txt
#SBATCH -e ../twas/logs/predict_expression/plink/compute_predict_expression_Amygdala_full_gene_1_%a_E.txt
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

module load conda_R/4.4

## List current modules
module list

# Define the full path to your single-column file list
FILELIST="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/Amygdala_gene_list_score_1.txt"

OUTPUT_DIR="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/expression_prediction"

# --- 1. Extract the score File Name ---

SCORE_FILENAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${FILELIST})

# Check if the line was empty (e.g., if array ID is past the end of the file)
if [ -z "${SCORE_FILENAME}" ]; then
    echo "Line ${SLURM_ARRAY_TASK_ID} is empty. Exiting."
    exit 0
fi

# --- 2. Create the Output Filename Pattern ---

SCORE_DIR="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/score/"

# Extract the gene name part (e.g., "gene_6479") by removing the ".wgt.RDat" suffix
OUTPUT_FILE_PREFIX=$(echo ${SCORE_FILENAME} | sed 's/_score\.txt$/_prediction/')

# Construct the full output path
OUTPUT_FILE="${OUTPUT_DIR}/${OUTPUT_FILE_PREFIX}"

echo "Processing input: ${SCORE_FILENAME}"
echo "Saving to output prefix to : ${OUTPUT_FILE}"

# --- 3. Run the Function ---

FULL_SCORE_PATH="${SCORE_DIR}${SCORE_FILENAME}"

echo "Attempting to use score file at : ${FULL_SCORE_PATH}" # <--- Add this for debugging


plink --bfile /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_unique_snps_bim/mdd_bpd_maf01.rsidAmygdala_uniqueSNPs --score ${FULL_SCORE_PATH} 1 2 4 --out ${OUTPUT_FILE}

echo "**** Job ends ****"
date

