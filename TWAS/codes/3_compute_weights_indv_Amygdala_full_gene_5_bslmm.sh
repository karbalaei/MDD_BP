#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=3_compute_weights_indv_Amygdala_full_gene_5_bslmm
#SBATCH -c 1
#SBATCH -o logs/Amygdala_gene/compute_weights_indv_Amygdala_full_gene_5_%a_O_bslmm.txt
#SBATCH -e logs/Amygdala_gene/compute_weights_indv_Amygdala_full_gene_5_%a_E_bslmm.txt
#SBATCH --mail-type=ALL
#SBATCH --array=1-4500%500

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
module load fusion_twas/github
module load conda_R


## List current modules
module list


# relative path for FILELIST
FILELIST=$(echo "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/input_ids_5.txt")


## File id and feature name
FEATURENUM=$(awk 'BEGIN {FS="\t"} {print $1}' ${FILELIST} | awk "NR==${SLURM_ARRAY_TASK_ID}")
FEATUREID=$(awk 'BEGIN {FS="\t"} {print $2}' ${FILELIST} | awk "NR==${SLURM_ARRAY_TASK_ID}")

## Define files
FILTBIM="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene/bim_files/Amygdala_gene_${FEATURENUM}/filtered_snps_Amygdala_gene_${FEATURENUM}"
TMPFILES="tmp_files_bslmm/gene_${FEATURENUM}"
OUTFILES="out_files_bslmm/gene_${FEATURENUM}"


cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/Amygdala_gene

Rscript ./FUSION.compute_weights_LIBD.R \
    --bfile ${FILTBIM} \
    --tmp ${TMPFILES} \
    --out ${OUTFILES} \
    --PATH_gemma /jhpce/shared/libd/core/fusion_twas/github/fusion_twas-master/gemma-0.98.5-linux-static-AMD64 \
    --model top1,blup,lasso,enet,bslmm   --hsq_p 1.0001 --verbose 1 --save_hsq


echo "**** Job ends ****"
date