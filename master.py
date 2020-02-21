#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-05 08:20:56


import os
import sys

from clump_snps import clump_snps
from match_snps import match_snps
from get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
from ld_expand_control_snp_set import ld_expand_control_snps
from combine_control_sets import combine_control_sets
# -----------
# MAIN
# -----------

analysis_name = 'test_ptb'
outputdir= "output"


# 1) run clump_snps.py
inputfile="test_data/input_gwas_file.tsv"
OutObj = clump_snps(inputfile, outputdir, analysis_name)


# 2) run snpnsap
match_OutObj =match_snps(OutObj.get('lead_snps_file'), 2, 'friends_ld05', outputdir, analysis_name)


# 3) run ld plink
csnps_file = match_OutObj.get('matched_snps_file')
ldexp_OutObj = get_ldsnps_for_control_snps(csnps_file, outputdir, analysis_name)



# 4) run ld_expand_control_snps()
matched_snps_file = match_OutObj.get('matched_snps_file')
input_snps_num_ld_file=OutObj.get("bin_by_ld_file")
ld_snps_for_ctrl_snps_dir=ldexp_OutObj.get('ld_snps_for_control_snps_dir')
set_index = 1
ld_match_Obj = ld_expand_control_snps(matched_snps_file, input_snps_num_ld_file, ld_snps_for_ctrl_snps_dir, outputdir, set_index, analysis_name)


# 5) combine all the control sets into a file
controlsets_dir=ld_match_Obj.root_dir
clumped_dir=OutObj.root_dir
combine_control_sets(controlsets_dir, clumped_dir)

