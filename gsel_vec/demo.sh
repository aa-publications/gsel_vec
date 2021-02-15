#!/bin/bash


outputpath=$(pwd)

input_gwas="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/gsel_vec/demo/input_data/bmi_40k.txt"
# input_gwas="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/gsel_vec/demo/input_data/bmi_small.test"


# ipython --pdb -c "%run  run_evo_sig.py bmi40k ${input_gwas} ${outputpath}"
ipython --pdb -c "%run  run_evo_sig_test.py bmi_small ${input_gwas} ${outputpath}"