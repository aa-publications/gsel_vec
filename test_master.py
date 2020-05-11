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

sys.path.append('/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec')
from clump_snps import clump_snps, clump_snp_list
from match_snps import match_snps
from get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
# from ld_expand_control_snp_set import ld_expand_control_snps
# from dev_ld_expand_all_control_sets import ld_expand_control_snps
from combine_control_sets import combine_control_sets
from dev_ld_expand_all_control_sets import ld_expand_all_control_snps
from check_ld_expanded_control_sets import check_ld_expanded_sets
from intersect_annotation import intersect_all_annotations
from helper_general import start_logger, safe_mkdir
import logging
from datetime import datetime


master_start = time.time()




###
###    dependecies
###

analysis_name ="test"
snpsnap_db_file="/dors/capra_lab/projects/gwas_allele_age_evolution/data/snpsnap_database/raw_download/eur/ld0.9_collection.tab.gz"
thous_gen_file = '/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/1kg/EUR.chr{}.phase3.nodups'
anno_dir="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/create_annotations/anno_dict"

anno_path_dict = {'argweave':      os.path.join(anno_dir, 'argweave_snpsnap_eur_ld0.1_collection.pickle'),
                'betascore':       os.path.join(anno_dir, 'betascore_snpsnap_eur_ld0.1_collection.pickle'),
                'fst_eas_afr':     os.path.join(anno_dir, 'fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle'),
                'fst_eur_afr':     os.path.join(anno_dir, 'fst_eur_afr_snpsnap_eur_ld0.1_collection.pickle'),
                'fst_eur_eas':     os.path.join(anno_dir, 'fst_eur_eas_snpsnap_eur_ld0.1_collection.pickle'),
                'gerp':            os.path.join(anno_dir, 'gerp_snpsnap_eur_ld0.1_collection.pickle'),
                'iES_Sabeti':      os.path.join(anno_dir, 'iES_Sabeti_snpsnap_eur_ld0.1_collection.pickle'),
                'linsigh':         os.path.join(anno_dir, 'linsigh_snpsnap_eur_ld0.1_collection.pickle'),
                'phastCon100':     os.path.join(anno_dir, 'phastCon100_snpsnap_eur_ld0.1_collection.pickle'),
                'phyloP100':       os.path.join(anno_dir, 'phyloP100_snpsnap_eur_ld0.1_collection.pickle'),
                'xpehh_afr2_eas':  os.path.join(anno_dir, 'xpehh_afr2_eas_snpsnap_eur_ld0.1_collection.pickle'),
                'xpehh_afr2_eur':  os.path.join(anno_dir, 'xpehh_afr2_eur_snpsnap_eur_ld0.1_collection.pickle'),
                'xpehh_eas_eur':   os.path.join(anno_dir, 'xpehh_eas_eur_snpsnap_eur_ld0.1_collection.pickle')}

two_sided_bool_dict = {'argweave': True,
             'betascore': False,
             'fst_eas_afr': False,
             'fst_eur_afr': False,
             'fst_eur_eas': False,
             'gerp': False,
             'iES_Sabeti': False,
             'linsigh': False,
             'phastCon100': False,
             'phyloP100': True,
             'xpehh_afr2_eas': False,
             'xpehh_afr2_eur': False,
             'xpehh_eas_eur': False}



###
###    inputs and outputs
###

gwas_summary_file="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/test/input_data/bmi_small.test" # gwas file
outputdir=f"/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/dev_{analysis_name}"
safe_mkdir(outputdir)


###
###    selection pipeline paramters
###

num_control_sets = 5
lead_snp_min_gwas_pvalue=0.00000005
ld_snps_min_gwas_pvalue=0.00000005
min_r2_to_clump=0.1
min_kb_from_index_to_clump=500
ld_expand_lead_snp_min_r2=0.9

ldbuds_r2_threshold="friends_ld09"
control_snps_ld_expand_r2=0.9

ld_thresholds = ['ld<=1.0']
summary_type='max'

###
###    start logger
###

logfile=os.path.join(outputdir,f"{analysis_name}.log")
logger = start_logger(logfile)

# -----------
#  START PIPELINE
# -----------

# clump snps based on LD
OutObj = clump_snps(gwas_summary_file, outputdir, thous_gen_file, lead_snp_min_gwas_pvalue=lead_snp_min_gwas_pvalue, ld_snps_min_gwas_pvalue=ld_snps_min_gwas_pvalue,
                    min_r2_to_clump=min_r2_to_clump, min_kb_from_index_to_clump=min_kb_from_index_to_clump, ld_expand_lead_snp_min_r2=ld_expand_lead_snp_min_r2)

# match snps
match_OutObj =match_snps(OutObj.get('lead_snp_ld_pairs_r2'), num_control_sets, ldbuds_r2_threshold, snpsnap_db_file, outputdir)

# get ld snps for control snps
csnps_file = match_OutObj.get('matched_snps_file')
ldexp_OutObj = get_ldsnps_for_control_snps(csnps_file, thous_gen_file, outputdir, control_snps_ld_expand_r2)

# ld expand control snps for each input/lead snp
lead_snps_ld_counts_file = OutObj.get("bin_by_ld_file")
gwas_snps_r2_file = OutObj.get("lead_snp_ld_pairs_r2")
matched_file = match_OutObj.get('matched_snps_file')
control_ld_dir = ldexp_OutObj.get('ld_snps_for_control_snps_dir')

ldexp_match_OutObj = ld_expand_all_control_snps(lead_snps_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, outputdir, ld_thresholds = ld_thresholds)


# check ld expanded control sets
ld_expanded_control_sets_file= ldexp_match_OutObj.get('ld_expanded_output')
ld_expanded_control_sets_r2_file= ldexp_match_OutObj.get('ld_r2_expanded_output')
combined_summary_df, mean_sd_ldscore_df = check_ld_expanded_sets(snpsnap_db_file, ld_expanded_control_sets_file , lead_snps_ld_counts_file, ld_expanded_control_sets_r2_file, ldbuds_r2_threshold, outputdir)


# intersect annotation
intersect_ouputs = intersect_all_annotations(anno_path_dict, two_sided_bool_dict, summary_type, ld_expanded_control_sets_file, outputdir)


logger.info("Done with master.py. Took {:.2f} minutes.".format( (time.time()-master_start)/60))

