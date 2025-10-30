#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=1G
#SBATCH --job-name=compute_score_file_Amygdala_full_gene_4
#SBATCH -c 1
#SBATCH -o ../twas/logs/predict_expression/score_file/compute_score_file_Amygdala_full_gene_4_%a_O.txt
#SBATCH -e ../twas/logs/predict_expression/score_file/compute_score_file_Amygdala_full_gene_4_%a_E.txt
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
module load conda_R/4.4

## List current modules
module list

# Define the full path to your single-column file list
FILELIST="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/Amygdala_exp_pred_gene_list_4.txt"

OUTPUT_DIR="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/score"

# --- 1. Extract the RDat File Name ---
# Since the input file has only one column (the RDat file name), we only need one awk command.
# The selected line contains the pattern "gene_6479.wgt.RDat"

RDAT_FILENAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${FILELIST})

# Check if the line was empty (e.g., if array ID is past the end of the file)
if [ -z "${RDAT_FILENAME}" ]; then
    echo "Line ${SLURM_ARRAY_TASK_ID} is empty. Exiting."
    exit 0
fi

# --- 2. Create the Output Filename Pattern ---
# The goal is to convert "gene_i.wgt.RDat" to "gene_i_score.txt"

RDAT_DIR="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/out_files/"

# Extract the gene name part (e.g., "gene_6479") by removing the ".wgt.RDat" suffix
GENE_PREFIX=$(echo ${RDAT_FILENAME} | sed 's/\.wgt\.RDat$//')

# Construct the full output path
OUTPUT_FILE="${OUTPUT_DIR}/${GENE_PREFIX}_score.txt"

echo "Processing input: ${RDAT_FILENAME}"
echo "Saving to output: ${OUTPUT_FILE}"

# --- 3. Run the Function ---
# Rscript make_score.R input > output

FULL_RDAT_PATH="${RDAT_DIR}${RDAT_FILENAME}"

echo "Attempting to use RDat file at: ${FULL_RDAT_PATH}" # <--- Add this for debugging


Rscript /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/make_score.R ${FULL_RDAT_PATH} > ${OUTPUT_FILE}

echo "**** Job ends ****"
date
