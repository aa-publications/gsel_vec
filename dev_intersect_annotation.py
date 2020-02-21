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

from functools import partial
from multiprocessing import Pool, cpu_count

from helper_general import Outputs

import logging
logger = logging.getLogger('main.{}'.format(__name__))

from datetime import datetime
DATE = datetime.now().strftime('%Y-%m-%d')



# -----------
# FUNCTIONS
# -----------
def pooled_overlap(long_df):
    # annotation overlap: pool all snps across control sets and summarize by lead snp
    anno_overlap_df  = long_df.groupby(['lead_snp'])['anno','is_na'].count().reset_index()
    anno_overlap_df['prop_coverage']= ((anno_overlap_df['anno']/anno_overlap_df['is_na'])*100).round(2)
    anno_overlap_df.sort_values('prop_coverage', inplace=True)

    return anno_overlap_df

def per_set_overlap(long_df):

    # annotation overlap: summarize by lead snp per control set then take mean
    by_control_set_overlap_df = long_df.groupby(['lead_snp','set'])['anno', 'is_na'].count().reset_index()
    by_control_set_overlap_df['prop_coverage'] = ((by_control_set_overlap_df['anno']/by_control_set_overlap_df['is_na'])*100).round(2)


    # exclude gwas snps
    by_control_summary_df = by_control_set_overlap_df.loc[by_control_set_overlap_df['set']!='ld_snp'].groupby(['lead_snp']).agg({'prop_coverage':['mean','std','size']}).reset_index()
    by_control_summary_df.columns = by_control_summary_df.columns.droplevel()
    by_control_summary_df.columns = ['lead_snp','mean_prop','std_prop','n_sets']
    by_control_summary_df.sort_values('mean_prop', inplace=True)

    return by_control_set_overlap_df

def groupby_mean(control_long_df, gwas_long_df, by):

    if by =='mean':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).mean().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).mean().unstack().reset_index()
    elif by =='median':
        gwas_df = gwas_long_df.groupby(['lead_snp', 'set']).median().unstack().reset_index()
        cntrl_df = control_long_df.groupby(['lead_snp', 'set']).median().unstack().reset_index()


    gwas_df.columns = gwas_df.columns.droplevel()
    gwas_df.columns.name=''
    gwas_df.columns = ['lead_snp', 'lead_snp_anno']

    cntrl_df.columns = cntrl_df.columns.droplevel()
    cntrl_df.columns.name=''
    cntrl_df.columns = ['lead_snp'] + cntrl_df.columns[1:].values.tolist()



    summary_df = pd.merge(gwas_df, cntrl_df, on='lead_snp', how='outer')

    return summary_df

def get_mean_median_by_lead_snp(long_df):

    control_long_df = long_df.loc[long_df['set'] != 'ld_snp', ['lead_snp', 'set','anno']].copy()
    gwas_long_df = long_df.loc[long_df['set'] == 'ld_snp', ['lead_snp', 'set','anno']].copy()


    # mean and meadian skip over np.nan
    # if nan is present in output, that means it tried to summarize over only np.nans
    mean_df = groupby_mean(control_long_df, gwas_long_df, 'mean')
    median_df = groupby_mean(control_long_df, gwas_long_df, 'mean')

    return mean_df, median_df

def intersect_annotation(anno_label_path_pair, matched_file, output_dir=None, anno_label=None):

    anno_label = anno_label_path_pair[0]
    anno_file = anno_label_path_pair[1]

    # load annotation and matched snps
    anno_dict = pickle.load(open(anno_file, 'rb'))
    match_df = pd.read_csv(matched_file, sep="\t")
    print(f"{anno_label} has {len(anno_dict):,} values.")

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

    # summarize by lead_snp
    mean_df, median_df = get_mean_median_by_lead_snp(long_df)

    #
    # # write
    # pooled_overlap_df.to_csv(pooled_overlap_output, sep="\t", index=False)
    # by_control_set_overlap_df.to_csv(by_control_set_overlap_output, sep="\t", index=False)
    # mean_df.to_csv(mean_output, sep="\t", index=False)
    # median_df.to_csv(median_output, sep="\t", index=False)
    # print("wrote intersection overlap, mean, median to {}".format(output_dir))



    return {'anno_label':anno_label, 'mean_df':mean_df, 'median_df':median_df, 'pooled_overlap': pooled_overlap_df, 'by_control_set_overlap':by_control_set_overlap_df}

def set_up_outputs(OutputObj, anno_label_list):


    for anno_label in anno_label_list:

        # for each type of annotation, have one directory
        OutputObj.add_output(anno_label, anno_label, mkdir=True, add_root=True)

        anno_root = os.path.join(OutputObj.root_dir, anno_label)

        # within that dir,  set up ouput files
        OutputObj.add_output('{}_pooled_overlap_output'.format(anno_label),
                                'pooled_{}_overlap.tsv'.format(anno_label) , custom_root=anno_root)
        OutputObj.add_output('{}_by_control_set_overlap_output'.format(anno_label),
                                'by_control_set_{}_overlap.tsv'.format(anno_label), custom_root=anno_root)
        OutputObj.add_output('{}_mean_output'.format(anno_label),
                                'mean_{}_per_lead_snp.tsv'.format(anno_label), custom_root=anno_root)
        OutputObj.add_output('{}_median_output'.format(anno_label),
                                'median_{}_per_lead_snp.tsv'.format(anno_label), custom_root=anno_root)


    return OutputObj

def intersect_all_annotations(anno_path_dict, matched_file, output_root, analysis_name):
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
    analysis_name : str
        a prefix for this analysis; will be used as a prefix for hte output directory that is created

    Returns
    -------
    intersect_ouputs: a list with multiple dictionary elements
        each element is a dictionary  corresponding to one slection annotion
        each dictionary contrains the annotaiton label and various outputs (see dictionary keys)


    """



    intersect_start = time.time()


    logger.info(f"Intersection of {len(anno_path_dict)} has begun...")
    anno_label_list, anno_paths = zip(*anno_path_dict.items())


    # set up outputs
    output_dir = os.path.join(output_root, '{}_annotation_intersection'.format(analysis_name))
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj, anno_label_list)


    # one thread per annotation
    num_threads = cpu_count()
    mstart = time.time()
    logger.info("Using {:,} cores to intersect {:,} annotations.".format(num_threads-1, len(anno_path_dict)))
    pool = Pool(num_threads-1)
    partial_intersect_anno = partial(intersect_annotation, matched_file=matched_file)
    # create label: filepath dictionary pairs
    anno_label_paths_pairs = [(anno_label,anno_path) for anno_label, anno_path in anno_path_dict.items()]
    intersect_ouputs = pool.map(partial_intersect_anno, anno_label_paths_pairs)

    pool.close()
    pool.join()
    logger.info("Done intersection all annotations. Took {:.2f} minutes.".format( (time.time() - mstart)/60))


    # write
    for anno_result in  intersect_ouputs:

        anno_label = anno_result['anno_label']
        mean_df = anno_result['mean_df']
        median_df = anno_result['median_df']
        pooled_overlap_df = anno_result['pooled_overlap']
        by_control_set_overlap = anno_result['by_control_set_overlap']


        mean_df.to_csv(OutObj.get('{}_mean_output'.format(anno_label)), sep="\t", index=False)
        median_df.to_csv(OutObj.get('{}_median_output'.format(anno_label)), sep="\t", index=False)
        pooled_overlap_df.to_csv(OutObj.get('{}_pooled_overlap_output'.format(anno_label)) , sep="\t", index=False)
        by_control_set_overlap.to_csv(OutObj.get('{}_by_control_set_overlap_output'.format(anno_label)) , sep="\t", index=False)


    logger.info(f"Wrote annotaiton intersection to: {output_dir}")
    logger.info(f"Done intersecting with annotations. Took { (time.time()-intersect_start)/60 :.2f} minutes.")


    return intersect_ouputs