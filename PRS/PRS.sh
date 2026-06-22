#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=3_Run_PRS
#SBATCH -c 4
#SBATCH -o logs/PRS_O.txt
#SBATCH -e logs/PRS_E.txt
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date

## Ensure log directory exists
mkdir -p logs

## Load dependencies
module load conda_R/4.3

## List current modules
module list

## Define files and directories
PRSICE_EXEC="/jhpce/shared/libd/core/PRSice_linux"
R_SCRIPT="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/PRS/PRSice.R"
BASE_GWAS="MDD_GWAS_base.txt"
TARGET_DATA="mdd_bpd_maf01"            # Prefix for .bed, .bim, .fam files
PHENOTYPE_FILE="cohort_phenotypes.txt"  # Linking Sample IDs to MDD/BP/Control status

# Run PRSice-2 
# Note: --thread matches the #SBATCH -c allocation
Rscript $R_SCRIPT \
    --prsice $PRSICE_EXEC \
    --base $BASE_GWAS \
    --target $TARGET_DATA \
    --pheno $PHENOTYPE_FILE \
    --snp SNP \
    --chr CHR \
    --bp BP \
    --A1 A1 \
    --A2 A2 \
    --stat OR \
    --pvalue P \
    --binary-target T \
    --out PRS_MDD_Results \
    --thread 4

echo "**** Job ends ****"
date