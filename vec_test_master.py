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


# sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_pipeline")
sys.path.append("/scratch/abraha1/gsel_/gsel_pipeline_vec")
from clump_snps import clump_snps
from match_snps import match_snps
from get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
from ld_expand_control_snp_set import ld_expand_control_snps
from combine_control_sets import combine_control_sets
from ld_expand_all_control_sets import ld_expand_all_control_snps
from check_ld_expanded_control_sets import check_ld_expanded_sets
from intersect_annotation import intersect_all_annotations

import logging
from datetime import datetime



def start_logger(log_file):

    logger = logging.getLogger('main')
    logger.setLevel(logging.DEBUG)

    # format1 = logging.Formatter("%(levelname)s - %(asctime)s %(name)s line#:%(lineno)d --> %(message)s")
    format1 = logging.Formatter("%(levelname)s [%(asctime)s] - %(name)s line:%(lineno)d -> %(message)s","%x %H:%M:%S")
    format2 = logging.Formatter("%(message)s")


    fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(format1)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(format2)

    logger.addHandler(fh)
    logger.addHandler(ch)

    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    logger.info(f'Logging to {log_file} start on {dt_string}')

    return logger

# -----------
# MAIN
# -----------
master_start = time.time()


# -----------
# parameters
# -----------
snpsnap_db_file="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_pipeline_vec/snpsnap_database/ld0.1_collection.tab.pickle"
thous_gen_file = '/scratch/abraha1/gsel_/gsel_pipeline_vec/1kg/EUR.chr{}.phase3.nodups'

#
inputfile="test/input_data/bmi_small.test"
analysis_name = 'giant_bmi_small'
outputdir= "test/bmi_small_vec/"

# inputfile="test/input_data/ptb_input_gwas_file.tsv"
# analysis_name = 'ptb'
# outputdir= "test/ptb/"

logfile=os.path.join(outputdir,"{}.log".format(analysis_name))
num_control_sets = 50
logger =start_logger(logfile)



# -----------
# start
# -----------

# CLUMP
OutObj = clump_snps(inputfile, outputdir, analysis_name, thous_gen_file)
ldbuds_r2_threshold="friends_ld05"

# MATCH
match_OutObj =match_snps(OutObj.get('clumped_snps_file'), num_control_sets, ldbuds_r2_threshold, snpsnap_db_file, outputdir, analysis_name)

# GET LD SNPS FOR CONTROL SNPS
csnps_file = match_OutObj.get('matched_snps_file')
ldexp_OutObj = get_ldsnps_for_control_snps(csnps_file, thous_gen_file, outputdir, analysis_name)

# LD EXPAND MATCHED REGIONS
lead_ld_counts_file = OutObj.get("bin_by_ld_file")
gwas_snps_r2_file = OutObj.get("lead_snp_ld_pairs_r2")
matched_file = match_OutObj.get('matched_snps_file')
control_ld_dir = ldexp_OutObj.get('ld_snps_for_control_snps_dir')
ldexp_match_OutObj = ld_expand_all_control_snps(lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, outputdir, analysis_name)


# check
ld_expanded_control_sets_file= ldexp_match_OutObj.get('ld_expanded_output')
ld_expanded_control_sets_r2_file= ldexp_match_OutObj.get('ld_r2_expanded_output')
combined_summary_df, mean_sd_ldscore_df = check_ld_expanded_sets(snpsnap_db_file, ld_expanded_control_sets_file , lead_ld_counts_file, ld_expanded_control_sets_r2_file, ldbuds_r2_threshold, outputdir, analysis_name)
logger.info("Done with master.py. Took {:.2f} minutes.".format( (time.time()-master_start)/60))


# intersect
anno_dir="/scratch/abraha1/gsel_/gsel_pipeline_vec/create_annotations/anno_dicts"
anno_path_dict = {'fst': os.path.join(anno_dir,'fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle'), 'fst2': os.path.join(anno_dir,'fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle')}
intersect_ouputs = intersect_all_annotations(anno_path_dict, ld_expanded_control_sets_file, outputdir, analysis_name)

intersect_all_annotations(anno_path_dict, matched_file, output_root, analysis_name)