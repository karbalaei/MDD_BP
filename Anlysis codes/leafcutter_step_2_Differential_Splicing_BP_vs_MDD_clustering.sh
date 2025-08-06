#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=leafcutter_step_2_Differential_Splicing_BP_vs_MDD_clustering
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_2_Differential_Splicing_BP_vs_MDD_clustering.txt
#SBATCH -e ../logs/e_leafcutter_step_2_Differential_Splicing_BP_vs_MDD_clustering.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"



ml conda

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/leafcutter/covariates #bacause of handeling intermedaite files!

python /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/codes/leafcutter_cluster_regtools.py -j ../junc_files/leafcutter_DS_juncfiles_covariates.txt -m 50 -l 500000 


echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.0.0
## available from http://research.libd.org/slurmjobs/
