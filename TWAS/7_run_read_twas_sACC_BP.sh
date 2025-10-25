#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=7_read_twas_BP_sACC
#SBATCH -c 1
#SBATCH -o logs/7_read_twas_BP_sACC_O.txt
#SBATCH -e logs/7_read_twas_BP_sACC_E.txt
#SBATCH --mail-type=ALL


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \$USER"
echo "Job id: \$SLURM_JOB_ID"
echo "Job name: \$SLURM_JOB_NAME"
echo "Hostname: \$HOSTNAME"
echo "Task id: \$SLURM_ARRAY_TASK_ID"


## Load dependencies
module load conda_R

## List current modules
module list

Rscript 7_read_twas_BP.R -r "sACC"

echo "**** Job ends ****"
date


