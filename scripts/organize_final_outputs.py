#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-02-13 09:22:14

import logging
import os
import sys
import time
from datetime import datetime

import numpy as np
import pandas as pd

from .helper_general import Outputs

DATE = datetime.now().strftime("%Y-%m-%d")


logger = logging.getLogger("main.{}".format(__name__))


# -----------
# FUNCTIONS
# -----------

def set_up_outputs(OutObj, anno_label_list):


    for anno_label in anno_label_list:

        # within that dir,  set up ouput files
        OutObj.add_output(
            "{}_pval_zscore_per_annotation".format(anno_label),
            "{}_z_score_pval.tsv".format(anno_label),
            add_root=True,
        )

    return OutObj

def load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict):
    anno_label_list, anno_paths = zip(*anno_path_dict.items())
    anno_pval_df = pd.DataFrame()
    anno_zscore_df = pd.DataFrame()
    anno_cover_df = pd.DataFrame()
    for anno_label in anno_label_list:
        pval_file_label = "{}_pvalue".format(anno_label)
        zscore_file_label = "{}_zscore".format(anno_label)
        per_control_file_label = "{}_by_control_set_overlap_output".format(anno_label)


        pval_file = intersectAnnoOutputObj.get(pval_file_label)
        zscore_file = intersectAnnoOutputObj.get(zscore_file_label)
        control_set_anno_cover_file = intersectAnnoOutputObj.get(per_control_file_label)


        pval_df = pd.read_csv(pval_file, sep="\t")
        zscore_df = pd.read_csv(zscore_file, sep="\t")
        control_set_anno_cover_df = pd.read_csv(control_set_anno_cover_file, sep="\t")

        anno_label = anno_label.replace("_", " ")
        pval_df['annotation'] = anno_label
        zscore_df['annotation'] = anno_label
        control_set_anno_cover_df['annotation'] = anno_label

        anno_pval_df = anno_pval_df.append(pval_df)
        anno_zscore_df = anno_zscore_df.append(zscore_df)
        anno_cover_df = anno_cover_df.append(control_set_anno_cover_df)

    return anno_pval_df, anno_zscore_df, anno_cover_df


def organize_final_outputs(intersectAnnoOutputObj, anno_path_dict, match_quality_per_lead_snp_df, match_summary_by_params_df, ldbuds_r2_threshold,output_root):


    logger.info(f"Finalizing outputs...")

    ###
    ### set up outputs
    ###
    anno_label_list, anno_paths = zip(*anno_path_dict.items())
    output_dir = os.path.join(output_root, "pval_zscore_per_anno")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj,anno_label_list)


    ###
    ###    tabulate matching quality per lead snp
    ###
    match_quality_per_lead_snp_df['prop_match_include'] = match_quality_per_lead_snp_df['mean_prop'] > 0.95

    ###
    ###    tabulate annotation intersection quality
    ###
    param_df= match_summary_by_params_df.loc[:, ['outside_sd_snp_maf', 'outside_sd_gene_count', 'outside_sd_dist_nearest_gene','outside_sd_{}'.format(ldbuds_r2_threshold)]].copy()
    match_summary_by_params_df['param_match_include'] = ~param_df.any(1)


    ###
    ###    write final outputs
    ###


    anno_pval_df, anno_zscore_df, anno_cover_df = load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict )
    # per annotation, evaluate intersection quality & write one file per annoation


    for annotation in anno_label_list:


        this_anno_pval_df = anno_pval_df.loc[anno_pval_df['annotation'] == annotation, ['lead_snp','pvalue','test_type', 'reject_h0_benj_hoch',  'corrected_pval_benj_hoch', 'annotation']].copy()
        this_anno_zscore_df = anno_zscore_df.loc[anno_zscore_df['annotation'] == annotation, ['lead_snp', 'z_score', 'lead_snp_anno', 'mean_controls']].copy()
        this_anno_cover_df = anno_cover_df.loc[anno_cover_df['annotation'] == annotation, ['lead_snp','mean_prop']].copy()



        this_anno_cover_df['anno_coverage_include'] = this_anno_cover_df['mean_prop'] > 0.95

        # I assume the two merging dataframes have the same lead snps
        zscore_pval_df = pd.merge(this_anno_pval_df, this_anno_zscore_df, on='lead_snp', how='inner')
        temp_include_df = pd.merge(zscore_pval_df, this_anno_pval_df, on='lead_snp', how='inner')
        temp_include_df =  pd.merge(temp_include_df, this_anno_cover_df.loc[:, ['lead_snp','anno_coverage_include']], on='lead_snp', how='inner')
        temp_include_df =  pd.merge(temp_include_df, match_quality_per_lead_snp_df.loc[:, ['lead_snp','prop_match_include']], on='lead_snp', how='inner')
        zscore_pval_flags_df =  pd.merge(temp_include_df, match_summary_by_params_df.loc[:, ['lead_snp','param_match_include']], on='lead_snp', how='inner')


        zscore_pval_flags_df['final_qc_include'] = zscore_pval_flags_df.loc[:, ['anno_coverage_include','prop_match_include','param_match_include']].all(1)


        cols_to_write =['lead_snp','z_score', 'lead_snp_anno', 'mean_controls','pvalue','test_type','reject_h0_bonfer','corrected_pval_benj_hoch', 'annotation', 'final_qc_include']
        to_write_df = zscore_pval_flags_df.loc[:, cols_to_write].copy()


        keep_write_df = to_write_df.loc[to_write_df['final_qc_include']==True]
        exclude_write_df = to_write_df.loc[to_write_df['final_qc_include']==False]

        # write
        to_write_df.to_csv(OutObj.get("{}_pval_zscore_per_annotation".format(annotation)), sep="\t", index=False)

    return OutObj