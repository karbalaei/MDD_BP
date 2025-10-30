#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=DS_leafcutter_step_4
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_4.txt
#SBATCH -e ../logs/e_leafcutter_step_4.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00



wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz

mv gencode.v25.annotation.gtf.gz ../leafcutter/covariates_results/

chmod u+x gtf2leafcutter.pl
chmod u+x prepare_results.R


./gtf2leafcutter.pl -o ../leafcutter/covariates_results/gencode_hg38 ../leafcutter/covariates_results/gencode.v25.annotation.gtf.gz





