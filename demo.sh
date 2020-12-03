#!/bin/bash


outputpath=$(pwd)
input_gwas="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/test/input_data/bmi_small.test"

# input_gwas="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/test/input_data/formatted_Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt"

python run_evo_sig.py demov1 $input_gwas $outputpath