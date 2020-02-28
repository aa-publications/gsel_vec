#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:10:00
#SBATCH --mem=100G
#SBATCH --output=inter_%A_%a.out
#SBATCH --array=1-5



# chr_dir="/dors/capra_lab/data/evolutionary_conservation/phastcons/hg19_100way/by_chr"



chr_dir="/dors/capra_lab/data/evolutionary_conservation/phyloP/phyloP_100way_hg19/by_chr"
snpsnap_file="/dors/capra_lab/projects/gwas_allele_age_evolution/data/snpsnap_database/chr${SLURM_ARRAY_TASK_ID}_snpsnap_ld0.1_EUR.bed.gz"


output_dir="/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/phyloP/snpsnap_intersect"


# awk -v OFS="\t" -F ":" '{print "chr"$1,$2,$2+1}' ${snpsnap_dir}"/chr"${chrom_num}_snpsnapld0.1_EUR.txt |bedtools intersect -a stdin -b ${chr_dir}"/chr"${chrom_num}.bed.gz -wb > ${output_dir}/intersect_snpsnapld0.1_EUR_chr${chrom_num}.tsv

echo "bedtools intersect -a $snpsnap_file -b ${chr_dir}"/chr"${SLURM_ARRAY_TASK_ID}.bed.gz -wb > ${output_dir}/intersect_snpsnapld0.1_EUR_chr${SLURM_ARRAY_TASK_ID}.tsv"
bedtools intersect -a $snpsnap_file -b ${chr_dir}"/chr"${SLURM_ARRAY_TASK_ID}.bed.gz -wb > ${output_dir}/intersect_snpsnapld0.1_EUR_chr${SLURM_ARRAY_TASK_ID}.tsv


