#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=leafcutter_step_6_Amygdala
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_6_Amygdala.txt
#SBATCH -e ../logs/e_leafcutter_step_6_Amygdala.txt
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00


set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"



module load conda_R/4.4


## List current modules for reproducibility
module list

## Edit with your job command
Rscript leafcutter_step_6_Amygdala.R


echo "**** Job ends ****"
date

