#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-05 08:20:56


import os
import sys
import time
import numpy as np
import pandas as pd

from functools import partial
from multiprocessing import Pool, cpu_count

sys.path.append("/scratch/abraha1/gsel_/gsel_vec")
from clump_snps import clump_snps, clump_snp_list
from match_snps import match_snps
from get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
from ld_expand_control_snp_set import ld_expand_control_snps
from combine_control_sets import combine_control_sets
from ld_expand_all_control_sets import ld_expand_all_control_snps
from check_ld_expanded_control_sets import check_ld_expanded_sets
from intersect_annotation import intersect_all_annotations
from helper_general import start_logger, safe_mkdir
import logging
from datetime import datetime


# %%

# -----------
# SET UP
# -----------

master_start = time.time()
analysis_name ="test_input_snps"

###
###    dependecies
###

snpsnap_db_file="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_pipeline_vec/snpsnap_database/ld0.1_collection.tab.pickle"
thous_gen_file = '/scratch/abraha1/gsel_/gsel_vec/1kg/EUR.chr{}.phase3.nodups'


###
###    inputs and outputs
###

inputfile="test/input_data/bmi_small.test" # gwas file
snps_list_file="/scratch/abraha1/gsel_/gsel_vec/test/input_data/rand_snp_list.txt" # list of snps
outputdir=f"/scratch/abraha1/gsel_/gsel_vec/test/{analysis_name}_gsel_vec"
safe_mkdir(outputdir)


###
###    selection pipeline paramters
###

clump_r2_threshold=0.9 # clumping r2 threshold
num_control_sets = 5


###
###    start logger
###

logfile=os.path.join(outputdir,f"{analysis_name}.log")
logger = start_logger(logfile)



# %%
# -----------
#  START PIPELINE
# -----------

# clump gwas inputs
# OutObj = clump_snps(inputfile, outputdir, thous_gen_file)

# clump list of snps
OutObj = clump_snp_list(snps_list_file, outputdir, thous_gen_file, clump_r2_threshold=0.9)

# MATCH
ldbuds_r2_threshold="friends_ld05"
match_OutObj =match_snps(OutObj.get('lead_snp_ld_pairs_r2'), num_control_sets, ldbuds_r2_threshold, snpsnap_db_file, outputdir)

# GET LD SNPS FOR CONTROL SNPS
csnps_file = match_OutObj.get('matched_snps_file')
ldexp_OutObj = get_ldsnps_for_control_snps(csnps_file, thous_gen_file, outputdir, clump_r2_threshold)

# LD EXPAND MATCHED REGIONS
lead_snps_ld_counts_file = OutObj.get("bin_by_ld_file")
gwas_snps_r2_file = OutObj.get("lead_snp_ld_pairs_r2")
matched_file = match_OutObj.get('matched_snps_file')
control_ld_dir = ldexp_OutObj.get('ld_snps_for_control_snps_dir')


ldexp_match_OutObj = ld_expand_all_control_snps(lead_snps_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, outputdir)


# CHECK
ld_expanded_control_sets_file= ldexp_match_OutObj.get('ld_expanded_output')
ld_expanded_control_sets_r2_file= ldexp_match_OutObj.get('ld_r2_expanded_output')
combined_summary_df, mean_sd_ldscore_df = check_ld_expanded_sets(snpsnap_db_file, ld_expanded_control_sets_file , lead_snps_ld_counts_file, ld_expanded_control_sets_r2_file, ldbuds_r2_threshold, outputdir)

# output_root=outputdir



# INTERSECT
anno_dir="/scratch/abraha1/gsel_/gsel_vec/create_annotations/anno_dicts"
anno_path_dict = {'fst': os.path.join(anno_dir,'fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle'), 'fst2': os.path.join(anno_dir,'fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle')}
intersect_ouputs = intersect_all_annotations(anno_path_dict, ld_expanded_control_sets_file, outputdir)

# output_root = outputdir
# matched_file=ld_expanded_control_sets_file


logger.info("Done with master.py. Took {:.2f} minutes.".format( (time.time()-master_start)/60))