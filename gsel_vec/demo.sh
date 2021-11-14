#!/bin/bash


outputpath=$(pwd)

input_gwas="demo/input_data/bmi_small.test"
# input_gwas="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/gsel_vec/demo/input_data/bmi_small.test"


# ipython --pdb -c "%run  run_evo_sig_test.py demo_test ${input_gwas} ${outputpath}"
python run_evo_sig.py demo_test ${input_gwas} ${outputpath}