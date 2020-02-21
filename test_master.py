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

from functools import partial
from multiprocessing import Pool, cpu_count


sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_pipeline")
from clump_snps import clump_snps
from match_snps import match_snps
from get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
from ld_expand_control_snp_set import ld_expand_control_snps
from combine_control_sets import combine_control_sets


import numpy as np
import pandas as pd

import logging


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

    logger.info(f"Running: {' '.join(sys.argv)}\n")
    logger.info(f'Logging to {log_file}')

    return logger

# -----------
# MAIN
# -----------
master_start = time.time()

# inputfile="test/input_data/formatted_Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt"
inputfile="test/input_data/bmi_small.test"
analysis_name = 'giant_bmi'
outputdir= "test/bmi_small/"
logfile=os.path.join(outputdir,"{}.log".format(analysis_name))
num_control_sets = 7
logger =start_logger(logfile)



# -----------
# start
# -----------


# 1) run clump_snps.py
OutObj = clump_snps(inputfile, outputdir, analysis_name)

# 2) run snpnsap
match_OutObj =match_snps(OutObj.get('clumped_snps_file'), num_control_sets, 'friends_ld05', outputdir, analysis_name)

# 3) run ld plink
csnps_file = match_OutObj.get('matched_snps_file')
ldexp_OutObj = get_ldsnps_for_control_snps(csnps_file, outputdir, analysis_name)


# 4) run ld_expand_control_snps()
matched_snps_file = match_OutObj.get('matched_snps_file')
input_snps_num_ld_file=OutObj.get("bin_by_ld_file")
ld_snps_for_ctrl_snps_dir=ldexp_OutObj.get('ld_snps_for_control_snps_dir')

# matched_snps_file = os.path.join(outputdir, '{}_matched_snps'.format(analysis_name),"matched_snps.tsv")
# input_snps_num_ld_file=os.path.join(outputdir, '{}_clump'.format(analysis_name), "lead_gwas_snps_with_ldsnps_counts.txt")
# ld_snps_for_ctrl_snps_dir=os.path.join(os.path.join(outputdir, '{}_get_ldsnps_for_control_snps'.format(analysis_name), "ld_snps_for_control_snps"))


# --- start multithreading ---
num_threads = cpu_count()-10


partial_ldex_csnps = partial(ld_expand_control_snps, matched_snps_file,
                                                     input_snps_num_ld_file,
                                                     ld_snps_for_ctrl_snps_dir,
                                                     outputdir,
                                                     analysis_name)


pool = Pool(num_threads-1)
# set_indexes=np.arange(1,num_control_sets+1,1)
ldex_outputs = pool.map(partial_ldex_csnps, range(1,num_control_sets))

pool.close()
pool.join()



# 5) combine all the control sets into a file
controlsets_dir=ldex_outputs[0].root_dir
clumped_dir=OutObj.root_dir
combine_control_sets(controlsets_dir, clumped_dir, outputdir)
logger.info("Done with master.py. Took {:.2} minutes.".format( (time.time()-master_start)/60))