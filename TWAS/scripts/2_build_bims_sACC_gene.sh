#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=2_build_bims_sACC_gene
#SBATCH -c 1
#SBATCH -o /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/logs/2_build_bims_sACC_gene_O.txt
#SBATCH -e /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/logs/2_build_bims_sACC_gene_E.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/

## Load dependencies
module load plink/1.90b
module load fusion_twas/github
module load conda_R

## List current modules
module list

## Compute weights for the given region/feature pair
Rscript ./hydeGoes_scSeq_mdd/MDD_vs_BP/twas/2_build_bims.R -c 10 -r "sACC"

echo "**** Job ends ****"
date
