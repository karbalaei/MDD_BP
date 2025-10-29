#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=1_filter_snps_Amygdala_gene
#SBATCH -c 1
#SBATCH -o /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/logs/1_filter_snps_Amygdala_O.txt
#SBATCH -e /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/logs/1_filter_snps_Amygdala_E.txt
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load plink/1.90b

module load conda_R

## List current modules for reproducibility
module list

#mkdir ./hydeGoes_scSeq_mdd/MDD_vs_BP/twas/logs

## Edit with your job command EDIT EDIT
Rscript ./hydeGoes_scSeq_mdd/MDD_vs_BP/twas/1_filter_snps_gene.R -r "Amygdala"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
