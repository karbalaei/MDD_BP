#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=8_plot_twas
#SBATCH -c 1
#SBATCH -o logs/plot_twas_O.txt
#SBATCH -e logs/plot_twas_E.txt
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

Rscript 8_generate_twas_plots.R

echo "**** Job ends ****"
date
