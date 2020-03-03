#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'



import os
import sys
import glob
import time
import argparse
import logging

import numpy as np
import pandas as pd
from subprocess import Popen, PIPE

from pdb import set_trace

def get_list_of_ld_snps(clump_df_path):
    '''return per row, the lead snp with a snp in ld with it by parsing the plink --clump output'''

    cdf = pd.read_csv(clump_df_path, sep="\s+")

    store_pairs = list()
    for row in cdf.iterrows():
        lead_snp = row[1].SNP
        sp2_row = row[1].SP2

        if sp2_row == 'NONE':
            store_pairs.append([lead_snp, 'NONE'])
        else:
            snps = [x[:-3] for x in sp2_row.split(",")]  # remove '(1)' at the end of each snp

            for ld_snp in snps:
                store_pairs.append([lead_snp, ld_snp])

    return store_pairs


def create_r2_dict(r2_file):
    ''' creates a nested dicionary of this structure: {snp_a:{snp_b1:r2_1, snpb2:r2_2}} '''

    tstart = time.time()
    logger = logging.getLogger('clump_snps.py')
    logger.info("Creating r2 dictionary.")

    r2_dict = dict()
    ld_df = pd.read_csv(r2_file, sep="\t")
    uniq_snp_a = ld_df.SNP_A.unique()

    for counter, snp_a in enumerate(uniq_snp_a):

        print(f"On {counter} of {len(uniq_snp_a)}")
        all_snp_b = ld_df.loc[ld_df['SNP_A'] == snp_a, 'SNP_B'].values
        all_r2 = ld_df.loc[ld_df['SNP_A'] == snp_a, 'R2'].values

        r2_dict[snp_a] = dict(zip(all_snp_b, all_r2))

    logger.info(f"Done!")
    return r2_dict


def bin_ldsnp_per_leadsnp(og_lead_ld_r2_df):
    ''' for each lead snp, find each snp in LD with it; calculate its r2 and then bin the pair by LD'''

    # -------
    # SET UP LD BINS
    # -------
    lead_ld_r2_df = og_lead_ld_r2_df.copy()
    lead_ld_r2_df.reset_index(inplace=True, drop=True)
    r2_bins = np.arange(11)/10

    # need to add one so that lenght of r2_bins can
    # be represented by the index of np.bincount
    num_slots_for_digitize_results = len(r2_bins) + 1

    # -------
    # SET UP FINAL DATA FRAME TO RETURN
    # -------
    col_names = ["ld<={}".format(x) for x in r2_bins]
    col_names.append("ld>1")
    final_count_df = pd.DataFrame(columns=col_names)


    # -------
    # FOR EACH LEAD SNP, BIN LD PARTNERS
    # -------


    # seperate lead snps w/o any LD partners
    lead_ld_r2_df.R2 = lead_ld_r2_df.R2.astype('object')  # coerce to make NONE comparison
    snps_without_ldbuddies = lead_ld_r2_df.loc[lead_ld_r2_df['R2'] == 'NONE', 'lead_snp'].values
    lead_ld_r2_df.drop(lead_ld_r2_df.loc[lead_ld_r2_df['R2'] == 'NONE', 'lead_snp'].index, inplace=True)

    # check if there is at least one lead snp
    # with ld partners to sort into LD BINS
    if lead_ld_r2_df.shape[0] != 0:
        lead_ld_r2_df.R2 = pd.to_numeric(lead_ld_r2_df.R2)

        # sort ld parterners and count number in each bin
        count_df = lead_ld_r2_df.groupby('lead_snp').R2.apply(lambda x: np.bincount(
            np.digitize(x, r2_bins, right=True), minlength=num_slots_for_digitize_results))

        rsID_index_labels = count_df.index
        binned_count_df = pd.DataFrame([list(x) for x in count_df.values],  columns=col_names, index=rsID_index_labels)
        final_count_df = final_count_df.append(binned_count_df)

    # create dataframe with 0 values for
    # LD bins for snps without LD buddies
    temp_data = [[0]*final_count_df.shape[1]]*len(snps_without_ldbuddies)
    no_ld_buds = pd.DataFrame(temp_data, index=snps_without_ldbuddies, columns=col_names)

    # add snps without any ld buddies to the final dataframe
    final_count_df = final_count_df.append(no_ld_buds)



    return final_count_df


def get_r2(ld_pairs, r2_dict):
    ''' return a dataframe with each lead:ld snp pair and their r2 value'''

    store_r2 = [np.nan]*len(ld_pairs)

    for index, snps in enumerate(ld_pairs):

        snp_a, snp_b = snps

        if snp_a == 'NONE' or snp_b == 'NONE':
            store_r2[index] = [snp_a, snp_b, 'NONE']
        else:

            try:
                r2 = r2_dict[snp_a][snp_b]
            except KeyError:
                try:
                    r2 = r2_dict[snp_b][snp_a]
                except KeyError:
                    raise KeyError("Input snp pair r2 value was not found!")

            store_r2[index] = [snp_a, snp_b, r2]

    r2df = pd.DataFrame(
        store_r2,
        columns=['lead_snp', 'ld_snp_pair', 'r2'])

    return r2df


def lead_snp_rsID_to_coord(concatenated_clumped_output_df):
    ''' dictionary of lead snps with this form: {rsID:chr_coord} '''
    df = concatenated_clumped_output_df

    rsID_coord_dict = dict(zip(df.SNP, df.CHR_BP))

    return rsID_coord_dict


def get_r2_for_lead_ld_snps(og_lead_ld_df, r2_df):
    """ return a df with lead and and ld snp and R2. lead snps w/o any SNPs in r2_df is included as "NONE for ld snps adn R2 """
    lead_ld_df = og_lead_ld_df.copy()
    lead_ld_df.reset_index(inplace=True, drop=True)


    # seperate the SNPS w/o and LD SNPS
    no_ld_lead_df = lead_ld_df[lead_ld_df['ld_snp'] == "NONE"].copy()
    no_ld_lead_df['R2'] = "NONE"
    lead_ld_df.drop(no_ld_lead_df.index, inplace=True)

    # merge the r2 data for INPUT SNPS
    lead_ld_df['lead_ld_pair'] = lead_ld_df['lead_snp'] + "_" + lead_ld_df['ld_snp']
    temp_m1 = pd.merge(lead_ld_df, r2_df.loc[:,['snpA_B','R2']], left_on='lead_ld_pair', right_on="snpA_B", how='inner')
    temp_m2 = pd.merge(lead_ld_df, r2_df.loc[:,['snpB_A','R2']], left_on='lead_ld_pair', right_on="snpB_A", how='inner')
    merged_df = pd.concat([temp_m1.drop(['snpA_B'], axis=1),temp_m2.drop(['snpB_A'], axis=1)], axis=0)

    # round the R2 to save memeory...
    merged_df.R2 = merged_df.R2.round(2)

    # add any snps in lead_ld_df that were not present in the r2_df
    #           >>> NOTE: a lead snp may have snps in LD based on the clumping parameters, but may not have any snps in the r2_df because that is determined by the ld_expand parameter...
    snps_wo_r2 = [x for x in set(lead_ld_df.lead_snp.unique()) -  set(merged_df.lead_snp.unique())]
    snps_wo_r2_df = pd.DataFrame( {'lead_snp':snps_wo_r2, 'ld_snp':['NONE']*len(snps_wo_r2), 'R2':['NONE']*len(snps_wo_r2)})
    no_ld_lead_df = no_ld_lead_df.append(snps_wo_r2_df)


    full_df = pd.concat([merged_df.drop('lead_ld_pair', axis=1), no_ld_lead_df], axis=0)

    assert set(full_df.lead_snp.unique()) == set(og_lead_ld_df.lead_snp.unique()), 'number of lead snps started with is not equal to the number finished after merging with r2_df'

    return full_df


# def write_gwas_snps_by_chr(gwas_df, gwas_stats_by_chr_dir, gwas_stats_by_chr_snps_only_dir, chromosome_column, basepair_column, gwasfilename, pvalue_column):
def write_gwas_sumstats_by_chr(gwas_df, gwas_stats_by_chr_dir, gwas_stats_by_chr_snps_only_dir, chromosome_column, basepair_column, gwasfilename, pvalue_column):

    gwas_df['variantID'] = gwas_df.loc[:, chromosome_column].map(str) + ":" + gwas_df.loc[:, basepair_column].map(str)


    gwas_stats_by_chr_files_dict = dict()
    gwas_snps_by_chr_files_dict = dict()
    for this_chr in gwas_df.loc[:, chromosome_column].unique():

        if (int(this_chr) == 23):
            continue

        output_file = os.path.join(gwas_stats_by_chr_dir, f"chr{this_chr}_{gwasfilename}")
        snps_only_file = os.path.join(gwas_stats_by_chr_snps_only_dir, f"chr{this_chr}_{gwasfilename}")

        # write to SNP and P value to file
        gwas_df.loc[gwas_df[chromosome_column] == this_chr, ["variantID", pvalue_column]].to_csv(output_file, sep="\t", index=False, header=["SNP", "P"])
        gwas_stats_by_chr_files_dict[this_chr] = output_file

        # write to SNP(chr:BP) to file
        gwas_df.loc[gwas_df[chromosome_column] == this_chr, ["variantID"]].to_csv(snps_only_file, sep="\t", index=False, header=None)
        gwas_snps_by_chr_files_dict[this_chr] = snps_only_file
    return gwas_stats_by_chr_files_dict, gwas_snps_by_chr_files_dict

def write_snp_list_by_chr(snp_list,  snps_by_chr_snps_only_dir):
    """ split input snp list by chromosome into different text files"""

    snps_df = pd.DataFrame(snp_list, columns=['chr_pos'])
    snps_df['chr'],snps_df['pos'] = snps_df.chr_pos.str.split(":", 1).str
    snps_df.sort_values('chr',inplace=True)


    snps_by_chr_files_dict = dict()
    for this_chr in snps_df.loc[:, 'chr'].unique():

        snps_only_file = os.path.join(snps_by_chr_snps_only_dir, "chr{}_from_snplist.txt")
        snps_df.loc[snps_df['chr'] == this_chr, 'chr_pos'].to_csv(snps_only_file.format(this_chr), sep="\t", index=False, header=False)
        snps_by_chr_files_dict[this_chr] =snps_only_file.format(this_chr)

    return snps_by_chr_files_dict
