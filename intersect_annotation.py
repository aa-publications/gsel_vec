#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import pickle
import time
import argparse

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests
from functools import partial
from multiprocessing import Pool, cpu_count

from helper_general import Outputs

import logging
logger = logging.getLogger('main.{}'.format(__name__))

from datetime import datetime
DATE = datetime.now().strftime('%Y-%m-%d')


# -----------
# CONSTANTS
# -----------


# -----------
# FUNCTIONS
# -----------
def pooled_overlap(long_df):
    # annotation overlap: pool all snps across control sets and summarize by lead snp
    count_df  = long_df.groupby(['lead_snp']).size().reset_index().rename(columns={0:'num_snps'})
    missing_df = long_df.groupby(['lead_snp'])['is_na'].sum().reset_index().rename(columns={'is_na':'num_missing_anno'})

    anno_coverage_df = pd.merge(count_df, missing_df, on='lead_snp', how='outer')
    anno_coverage_df['prop_coverage']= 100-((anno_coverage_df['num_missing_anno']/anno_coverage_df['num_snps'])*100).round(1)

    anno_coverage_df.sort_values('prop_coverage', inplace=True)

    return anno_coverage_df

def per_set_overlap(long_df):

    # annotation overlap: summarize by lead snp per control set then take mean
    count_df  = long_df.loc[long_df['set']!='ld_snp'].groupby(['lead_snp', 'set']).size().reset_index().rename(columns={0:'num_snps'})
    missing_df = long_df.loc[long_df['set']!='ld_snp'].groupby(['lead_snp', 'set'])['is_na'].sum().reset_index().rename(columns={'is_na':'num_missing_anno'})

    # prop_coverage  ==  for each lead snp (ld expanded) and control set pair --> # of control snps with an annotation value
    anno_coverage_df = pd.merge(count_df, missing_df, on=['lead_snp', 'set'], how='outer')
    anno_coverage_df['prop_coverage']= 100-((anno_coverage_df['num_missing_anno']/anno_coverage_df['num_snps'])*100).round(1)


    # prop_coverage --> take the mean of prop_coverage across all the control sets
    # exclude gwas snps
    by_control_summary_df = anno_coverage_df.groupby(['lead_snp']).agg({'prop_coverage':['mean','std','count']}).reset_index()
    by_control_summary_df.columns = by_control_summary_df.columns.droplevel()
    by_control_summary_df.columns = ['lead_snp','mean_prop','std_prop','non_na_sets']
    by_control_summary_df.sort_values('mean_prop', inplace=True)

    return by_control_summary_df

def grouby_stat(control_long_df, gwas_long_df, by):

    if by =='mean':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).mean().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).mean().unstack().reset_index()
    elif by =='median':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).median().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).median().unstack().reset_index()
    elif by =='max':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).max().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).max().unstack().reset_index()
    elif by =='min':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).min().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).min().unstack().reset_index()


    gwas_df.columns = gwas_df.columns.droplevel()
    gwas_df.columns.name=''
    gwas_df.columns = ['lead_snp', 'lead_snp_anno']

    cntrl_df.columns = cntrl_df.columns.droplevel()
    cntrl_df.columns.name=''
    cntrl_df.columns = ['lead_snp'] + cntrl_df.columns[1:].values.tolist()



    summary_df = pd.merge(gwas_df, cntrl_df, on='lead_snp', how='outer')

    return summary_df

def get_summary_by_lead_snp(long_df, summary_by):

    control_long_df = long_df.loc[long_df['set'] != 'ld_snp', ['lead_snp', 'set','anno']].copy()
    gwas_long_df = long_df.loc[long_df['set'] == 'ld_snp', ['lead_snp', 'set','anno']].copy()

    assert summary_by in ['mean','median','max','min'], 'summary_by values is not valid'

    # mean and meadian skip over np.nan
    # if nan is present in output, that means it tried to summarize over only np.nans
    # mean_df = grouby_stat(control_long_df, gwas_long_df, 'mean')
    # median_df = grouby_stat(control_long_df, gwas_long_df, 'median')
    summary_by_df = grouby_stat(control_long_df, gwas_long_df, summary_by)

    return summary_by_df

def intersect_annotation(anno_label_path_pair, matched_file, two_sided_bool_dict, summary_type):

    anno_label = anno_label_path_pair[0]
    anno_file = anno_label_path_pair[1]

    # load annotation and matched snps
    anno_dict = pickle.load(open(anno_file, 'rb'))
    match_df = pd.read_csv(matched_file, sep="\t")
    logger.debug(f"{anno_label} has {len(anno_dict):,} values.")

    # convert to long df
    # note: 'ld_snp' refers to the gwas snp that is the lead or its ld snp
    control_cols = [col for col in  match_df.columns if col.startswith("Set_")]
    long_df = pd.melt(match_df, id_vars=['lead_snp'], value_vars=['ld_snp']+ control_cols, var_name='set', value_name='snp')

    long_df['anno'] = long_df.snp.map(anno_dict)
    long_df['is_na'] = long_df.anno.isnull()



    # annotation overlap: pool all snps across control sets and summarize by lead snp
    pooled_overlap_df = pooled_overlap(long_df)

    # annotation overlap: summarize by lead snp per control set then take mean
    by_control_set_overlap_df = per_set_overlap(long_df)


    #pval
    # summarize by lead_snp
    summary_by_df = get_summary_by_lead_snp(long_df, summary_type)
    summary_by_pval_df, summary_by_z_score_df = create_pval_zscore_df(summary_by_df, anno_label, two_sided_bool_dict)


    return {'anno_label':anno_label, 'summary_by_df':summary_by_df,
            'summary_by_pval_df':summary_by_pval_df, 'summary_by_zscore_df':summary_by_z_score_df,
            'summary_type':summary_type,
            'pooled_overlap': pooled_overlap_df, 'by_control_set_overlap':by_control_set_overlap_df}

def set_up_outputs(OutputObj, anno_label_list, summary_type):


    for anno_label in anno_label_list:

        # for each type of annotation, have one directory
        OutputObj.add_output(anno_label, anno_label, mkdir=True, add_root=True)

        anno_root = os.path.join(OutputObj.root_dir, anno_label)

        # within that dir,  set up ouput files
        OutputObj.add_output('{}_pooled_overlap_output'.format(anno_label),
                                'pooled_{}_overlap.tsv'.format(anno_label) , custom_root=anno_root)
        OutputObj.add_output('{}_by_control_set_overlap_output'.format(anno_label),
                                'by_control_set_{}_overlap.tsv'.format(anno_label), custom_root=anno_root)
        OutputObj.add_output('{}_{}_output'.format(anno_label, summary_type),
                                '{}_{}_per_lead_snp.tsv'.format(summary_type, anno_label), custom_root=anno_root)
        # OutputObj.add_output('{}_{}_output'.format(anno_label, summary_type),
        #                         '{}_{}_per_lead_snp.tsv'.format(summary_type, anno_label), custom_root=anno_root)
        OutputObj.add_output('{}_pvalue'.format(anno_label),
                                'pvalue_{}_per_lead_snp.tsv'.format(anno_label), custom_root=anno_root)
        OutputObj.add_output('{}_zscore'.format(anno_label),
                                'zscore_{}_per_lead_snp.tsv'.format(anno_label), custom_root=anno_root)

        OutputObj.add_output('{}_enrichment_summary'.format(anno_label),
                                '{}_enrichment_summary.tsv'.format(anno_label), custom_root=anno_root)


    return OutputObj

def create_pval_zscore_df(summary_df, anno_label, two_sided_bool_dict):

    # wide_df = long_df.pivot(index='lead_snp', columns='set_num', values='Set_')


    control_cols = [col for col in  summary_df.columns if col.startswith("Set_")]

    z_score_df = summary_df.copy()
    z_score_df['mean_controls'] = z_score_df.loc[:,control_cols].mean(1)
    z_score_df['std_controls'] = z_score_df.loc[:,control_cols].std(1)
    z_score_df['z_score'] = (z_score_df.lead_snp_anno - z_score_df.mean_controls)/z_score_df.std_controls
    z_score_df['n_controls'] = z_score_df.loc[:, control_cols].apply(lambda x: np.sum(~x.isnull()), axis=1)
    z_score_df = z_score_df.loc[:, ['lead_snp', 'lead_snp_anno', 'z_score', 'mean_controls', 'std_controls']].copy()
    z_score_df = z_score_df.round(3)


    two_tailed_bool=two_sided_bool_dict[anno_label]

    # calc p-value
    partial_pval = partial(calc_pval_and_summary, control_cols=control_cols, two_tailed_bool=two_tailed_bool)
    pval_df = summary_df.apply(partial_pval, axis=1, result_type='expand')
    pval_df.columns = ['lead_snp','pvalue', 'test_type', 'lead_snp_anno', 'input_percentile', 'num_control_snps']



    no_na_pval_df = pval_df[~pval_df.isnull().any(1)].copy()
    na_pval_lead_snps = pval_df.loc[pval_df.isnull().any(1), 'lead_snp'].values.tolist()

    pvalues = no_na_pval_df.pvalue

    # multiple testing correction: bonferroni
    # hypo_reject_bonfer, pval_corrected_bonfer, _, _ = multipletests(pvalues, alpha=0.05, method='bonferroni')


    # multiple testing correction: benjamini-hochberg
    hypo_reject_bh, pval_corrected_bh, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')

    # no_na_pval_df['reject_h0_bonfer'] = hypo_reject_bonfer
    # no_na_pval_df['corrected_pval_bonf'] = pval_corrected_bonfer
    no_na_pval_df['reject_h0_benj_hoch'] = hypo_reject_bh
    no_na_pval_df['corrected_pval_benj_hoch'] = pval_corrected_bh
    no_na_pval_df.sort_values('corrected_pval_benj_hoch', inplace=True, ascending=True)
    no_na_pval_df.reset_index(inplace=True, drop=True)




    na_df = pd.DataFrame({'lead_snp':na_pval_lead_snps, 'pvalue':np.nan, 'test_type':np.nan, 'lead_snp_anno':np.nan, 'input_percentile':np.nan, 'num_control_snps':np.nan,
                   'reject_h0_benj_hoch':np.nan,
                   'corrected_pval_benj_hoch':np.nan})


    all_pval_df = pd.concat((na_df, no_na_pval_df), axis=0)

    return all_pval_df, z_score_df

def calc_pval_and_summary(values, control_cols, two_tailed_bool):

    np.seterr(all='warn')


    # note will remove missing values or na before calculating p-value
    lead_snp = values.lead_snp
    input_snp_value = values.lead_snp_anno


    control_snp_values = values[control_cols]
    na_removed_control_values = control_snp_values[~control_snp_values.isnull()]

    # A percentileofscore of 80% means that 80% of values are less than or equal to the provided score.
    input_percentile = stats.percentileofscore(na_removed_control_values, input_snp_value,  kind='weak')




    # modify pval calc based on one tailed or two tailed test; have a default also
    if two_tailed_bool:

        pval = np.sum(na_removed_control_values >= input_snp_value)/len(na_removed_control_values)

        if pval > 0.5:
            pval = np.sum(na_removed_control_values <= input_snp_value)/len(na_removed_control_values)

        pval = pval*2
        test_type = 'two-tailed'
    else:

        pval = np.sum(na_removed_control_values >= input_snp_value)/len(na_removed_control_values)
        test_type = 'greater-than'


    # old way of doing on
    # if np.isnan(input_snp_value):
    #     pval=np.nan
    # else:
    #     pval = np.sum(na_removed_control_values > input_snp_value)/len(na_removed_control_values)



    num_control_snps = len(na_removed_control_values)



    summary = [lead_snp, pval, test_type, input_snp_value, input_percentile, num_control_snps]

    return summary

def intersect_all_annotations(anno_path_dict, two_sided_bool_dict, summary_type,  matched_file, output_root):
    """For matched control sets, intersect them with a selection annotation and return the mean and median for each locus.
            * each locus is defined by lead snp and LD snps (as defined earlier in the pipeline)

    Parameters
    ----------
    anno_path_dict : dict
        a dictionary key value pair:  'annotation_label' and 'annotation_file_path'. NOTE: path should be a pickled dictionary
    matched_file : str
        full path to the table with ld expanded lead snps and  control sets
    output_root : str
        path to which output directory for this analysis should be created


    Returns
    -------
    intersect_ouputs: a list with multiple dictionary elements
        each element is a dictionary  corresponding to one slection annotion
        each dictionary contrains the annotaiton label and various outputs (see dictionary keys)


    """



    intersect_start = time.time()


    logger.info(f"Intersecting {len(anno_path_dict)} selection measures.")
    anno_label_list, anno_paths = zip(*anno_path_dict.items())

    ###
    ### set up outputs
    ###

    output_dir = os.path.join(output_root, 'selection_intersected_matched_sets')
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj, anno_label_list, summary_type)


    ###
    ### one thread per annotation intersectiong
    ###

    num_threads = cpu_count()
    mstart = time.time()
    num_processes = 4
    logger.info("Using {:,} cores to intersect {:,} annotations.".format(num_processes, len(anno_path_dict)))
    # pool = Pool(processes=num_processes, maxtasksperchild=10)
    pool = Pool(processes=1, maxtasksperchild=1)
    partial_intersect_anno = partial(intersect_annotation, matched_file=matched_file, two_sided_bool_dict=two_sided_bool_dict, summary_type=summary_type)

    # create label: filepath dictionary pairs
    anno_label_paths_pairs = [(anno_label,anno_path) for anno_label, anno_path in anno_path_dict.items()]
    intersect_ouputs = pool.map(partial_intersect_anno, anno_label_paths_pairs)

    pool.close()
    pool.join()
    logger.info("Done intersection on all annotations. Took {:.2f} minutes.".format( (time.time() - mstart)/60))


    # unpack and write
    for anno_result in  intersect_ouputs:

        anno_label = anno_result['anno_label']

        summary_by_df = anno_result['summary_by_df']
        summary_by_pval_df = anno_result['summary_by_pval_df']
        summary_by_zscore_df = anno_result['summary_by_zscore_df']
        summary_type = anno_result['summary_type']
        pooled_overlap_df = anno_result['pooled_overlap']
        by_control_set_overlap = anno_result['by_control_set_overlap']

        # mean_pvalue_df = anno_result['mean_pval_df']
        # mean_zscore_df = anno_result['mean_zscore_df']
        # median_pvalue_df = anno_result['median_pval_df']
        # median_zscore_df = anno_result['median_zscore_df']



        combined_df = pd.merge( pd.merge(by_control_set_overlap.loc[:, ['lead_snp','mean_prop','std_prop']],
                    summary_by_df.loc[:, ['lead_snp','lead_snp_anno']], on='lead_snp', how='outer'),
                    summary_by_pval_df.loc[:, ['lead_snp','pvalue','reject_h0_benj_hoch','corrected_pval_benj_hoch']], on='lead_snp', how='outer')

        combined_df['summary_type'] = summary_type
        combined_df.mean_prop = combined_df.mean_prop.round(2)
        combined_df.std_prop = combined_df.std_prop.round(2)
        combined_df.lead_snp_anno = combined_df.lead_snp_anno.round(4)


        combined_df.to_csv(OutObj.get('{}_enrichment_summary'.format(anno_label)), sep="\t", index=False)
        summary_by_df.to_csv(OutObj.get('{}_{}_output'.format(anno_label, summary_type)), sep="\t", index=False)
        pooled_overlap_df.to_csv(OutObj.get('{}_pooled_overlap_output'.format(anno_label)) , sep="\t", index=False)
        by_control_set_overlap.to_csv(OutObj.get('{}_by_control_set_overlap_output'.format(anno_label)) , sep="\t", index=False)



        # add the type of summary stat used for selection measures ...
        summary_by_pval_df.to_csv(OutObj.get('{}_pvalue'.format(anno_label)) , sep="\t", index=False, na_rep="NaN")
        summary_by_pval_df['summary_type_over_input_region']=summary_type
        summary_by_zscore_df.to_csv(OutObj.get('{}_zscore'.format(anno_label)) , sep="\t", index=False, na_rep="NaN")
        summary_by_zscore_df['summary_type_over_input_region']=summary_type




    logger.info(f"Wrote annotaiton intersection to: {output_dir}")
    logger.info(f"Done intersecting with annotations. Took { (time.time()-intersect_start)/60 :.2f} minutes.")


    return intersect_ouputs
