#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=running_burden_test_Amygdala_BP
#SBATCH -c 1
#SBATCH -o logs/6_running_burden_test_Amygdala_BP_E.txt
#SBATCH -e logs/6_running_burden_test_Amygdala_BP_O.txt
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --mail-user karbalaei@jhmi.edu

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

#find ./Results/Amygdala/ -type f -name "*_isoTWAS.RDS" | sed 's#./Results/Amygdala/##g; s#_isoTWAS.RDS##g' > Amygdala_RDS_list

## List current modules
module list

Rscript 6_running_burden_test_BP.R -r "Amygdala"

echo "**** Job ends ****"
date