#!/bin/bash

## Usage:
# sh apply_weights.sh

mkdir -p logs

for region in Amygdala sACC
do

    for feature in gene
    # for feature in gene exon jxn tx
    do

        # set of summary stats
        for summstats in bip
        do

            SHORT="6_apply_weights_full_BP_${region}_${feature}_${summstats}"

            # Construct shell file
            echo "Creating script for region ${region} at the ${feature} level for ${summstats}"

            # Use 'cat <<-EOF' to strip leading tabs and newlines.
            cat > .${SHORT}.sh <<-EOF
#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=run_process_gwas
#SBATCH -c 1
#SBATCH -o logs/${SHORT}_O.txt
#SBATCH -e logs/${SHORT}_E.txt
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
module load fusion_twas/github
module load conda_R

## List current modules
module list

## Choose the correct GWAS summary statistics file
if [ "${summstats}" == "bip" ]
then
    summstatsfile="/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/hg38/bip2024_multianc_no23andMe_CLEAN_hg38.txt"
else
    echo "Unexpected ${summstats} input"
fi

## Apply weights for the given region/feature pair and the given GWAS summary statistics
mkdir -p "${region}_${feature}/${summstats}"


for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome \${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/libd/core/fusion_twas/github/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats \${summstatsfile} \
    --weights /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/${region}_${feature}/${region}_${feature}.pos \
    --weights_dir /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/${region}_${feature}/ \
    --ref_ld_chr /dcs04/lieber/lcolladotor/annotationFiles_LIBD001/fusion_twas_LDREF_hg38/1000G.EUR. \
    --chr \${chr} \
    --out ${region}_${feature}/${summstats}/${summstats}.\${chr}.dat

    echo ""
    echo "making plots for chromosome \${chr}"
    date
    echo ""

## companion post-processing step (plots only for genes)

if [ "\$feature" == "gene" ]
then  
    Rscript /jhpce/shared/libd/core/fusion_twas/github/fusion_twas-master/FUSION.post_process.R \
        --sumstats /\${summstatsfile} \
        --input ${region}_${feature}/${summstats}/${summstats}.\${chr}.dat \
        --out ${region}_${feature}/${summstats}/${summstats}.\${chr}.analysis \
        --ref_ld_chr /dcs04/lieber/lcolladotor/annotationFiles_LIBD001/fusion_twas_LDREF_hg38/1000G.EUR. \
        --chr \${chr} \
        --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr \
        --glist_path "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/glist-hg38"
else
    Rscript /jhpce/shared/libd/core/fusion_twas/github/fusion_twas-master/FUSION.post_process.R \
        --sumstats /\${summstatsfile} \
        --input ${region}_${feature}/${summstats}/${summstats}.\${chr}.dat \
        --out ${region}_${feature}/${summstats}/${summstats}.\${chr}.analysis \
        --ref_ld_chr /dcs04/lieber/lcolladotor/annotationFiles_LIBD001/fusion_twas_LDREF_hg38/1000G.EUR. \
        --chr \${chr} \
        --locus_win 100000 --verbose 2 --plot_corr
fi

done

echo "**** Job ends ****"
date
EOF
            call="sbatch .${SHORT}.sh"
            echo $call
            $call
        done
    done
done

