#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=7_FDR_calculaton_isotwas
#SBATCH -c 1
#SBATCH -o logs/7_FDR_calculaton_isotwas_%a_O.txt
#SBATCH -e logs/7_FDR_calculaton_isotwas_%a_E.txt
#SBATCH --mail-type=ALL
#SBATCH --array=1-4
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

#echo -e "Amygdala\tBP\nsACC\tBP\nAmygdala\tMDD\nsACC\tMDD"  > jobs
#echo -e "Amygdala\tBP\nsACC\tBP\nAmygdala\tMDD\nsACC\tMDD" | column -t -s $'\t' > jobs

# relative path for FILELIST
#job_list="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/jobs"

# Define the data as arrays
regions=("Amygdala" "sACC" "Amygdala" "sACC")
diagnoses=("BP" "BP" "MDD" "MDD")

# The output file
job_list="jobs"

# Create the jobs file using a for loop
rm -f "${job_list}" # Remove the file if it exists

for i in "${!regions[@]}"; do
  echo -e "${regions[$i]}\t${diagnoses[$i]}" >> "${job_list}"
done

# Optional: Verify the content of the jobs file
cat "${job_list}"

## File id and feature name

region=$(awk '{print $1}' "${job_list}" | awk "NR==${SLURM_ARRAY_TASK_ID}")
diag=$(awk '{print $2}' "${job_list}" | awk "NR==${SLURM_ARRAY_TASK_ID}")

#debugging command, check what the names are.
echo "region: $region"
echo "diagnosis: $diag"

## Edit with your job command EDIT EDIT
Rscript 7_FDR_calculaton_isotwas.R \
	--region "${region}"\
	--diag "${diag}" 

echo "**** Job ends ****"
date

