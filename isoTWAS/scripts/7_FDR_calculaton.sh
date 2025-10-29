#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=7_FDR_calculaton_isotwas
#SBATCH -c 1
#SBATCH -o logs/7_FDR_calculaton_isotwas_%a_O.txt
#SBATCH -e logs/7_FDR_calculaton_isotwas_%a_E.txt
#SBATCH --mail-type=ALL
#SBATCH --array=1-8
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
module load conda_R

## List current modules for reproducibility
module list

## File id and feature name

input=$(find ./Results/ -name "*.txt" | sed 's#./Results/##g' | awk "NR==${SLURM_ARRAY_TASK_ID}")
        
output=$(find ./Results/ -name "*.txt" | sed -e 's/BurdenTest/Final/g' -e 's/\.txt$/.xlsx/g' -e 's#./Results/##g' | awk "NR==${SLURM_ARRAY_TASK_ID}")

#debugging command, check what the names are.
echo "input file: $input"
echo "output file: $output"

## Edit with your job command EDIT EDIT
Rscript 7_FDR_calculaton.R \
	--file "${input}"\
	--out "${output}" 

echo "**** Job ends ****"
date

