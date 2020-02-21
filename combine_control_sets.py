#!/bin/python
# This script will combined all the control snps with ld snps into one text file. One column per control set.
#       INPUTS:
#           * see argparse for details
#
#       OUTPUTS:
#           * returns dataframe with one row per snp
#           * first 3 columns are 1) lead gwas snps defining locus 2) gwas input snp 3) lead_snp column is True is that row is a lead snp
#           * remaining columns are control sets
#
#       Note:
#           * where locus_lead_snp == gwas_snp means these are lead snps
#           * final control set written removes control sets that are missing > 5% of required snps
#
#
#       DEPENDENCIES:
#           * need to specify what r2 threshold to consider for number of control snps in ld with lead snp
#
# Abin Abraham
# created on: 2018-09-16 14:24:45

import os
import sys
import glob
import time
import argparse

import numpy as np
import pandas as pd



import pdb

import logging
logger = logging.getLogger('main.{}'.format(__name__))


# -----------
# FUNCTIONS
# -----------

def parse_input_args():
    parser = argparse.ArgumentParser(description='Combines control snps with ld snps into one file.')

    parser.add_argument('control_sets_dir', action='store',
                        type=str,
                        help='path to output from expand_control_set.py')

    parser.add_argument('clumped_dir', action='store',
                        type=str,
                        help='path to output from clump_snps.py')

    args = parser.parse_args()
    controlsets_dir = args.control_sets_dir
    clumped_dir = args.clumped_dir


    return control_sets_dir, clumped_dir


def add_control_column(set_file, input_df):
    '''Lines up control lead and ld snps with input lead and ld snps.

        inputs:
            * set_file
                - full path to control_set
            * input_df
                - should have columns: 'locus_lead_snp','input_snp','lead_snp' (True/False) at minimum

        output:
            * input_df
                - will append control set column to input_df
            * store_mismatches
                - any instance where more than required ld control snps exists

        details:
            * cells with np.nan in output means not enough control ld snps were found
            * errors will be saved in store_mismatches

    '''

    set_name = os.path.split(set_file)[1].split("-")[0]

    # create gwas_snp:control_snp dictionary from **control set data**
    control_df = pd.read_csv(set_file, sep="\t")
    lead_pairs = dict(zip(control_df.loc[control_df['lead_snp'] == True, 'gwas_snp'],
                          control_df.loc[control_df['lead_snp'] == True, 'control_snp']))

    store_mismatches = []
    for gwas_snp, control_snp in lead_pairs.items():

        # print("Lining up lead and ld control snps with input gwas snps.\ngwas_lead = {}\t control_lead = {}\n".format(gwas_snp, control_snp))

        # add control lead snps to lead input snp rows
        input_df.loc[(input_df['locus_lead_snp'] == gwas_snp) & (
            input_df['input_snp'] == gwas_snp), set_name] = control_snp

        # add control ld snps to remaining rows
        ld_controls = control_df.loc[(control_df['gwas_snp'] == gwas_snp) & (
            control_df['lead_snp'] == False)].control_snp.values

        open_slots = input_df.loc[(input_df['locus_lead_snp'] == gwas_snp) & (
            input_df['lead_snp'] == False) & (input_df[set_name].isnull())]

        num_control_snps = len(ld_controls)
        num_slots = open_slots.shape[0]

        if num_control_snps > num_slots:
            # don't fills any of them
            store_mismatches.append('Num control ld snps > num LD snps based on gwas input.\
                    FOR: {}, gwas_snp {}, control_snp {}'.format(set_name, gwas_snp, control_snp))

        else:

            indicies_to_fill = open_slots.index[0:num_control_snps]
            input_df.loc[indicies_to_fill, set_name] = ld_controls

    return input_df, store_mismatches


def create_input_df(clumped_file, input_r2_file, upperinclusivebound, lowerexclusivebound):
    '''Creates table with lead snp paired with ld snps and filters based on r2 threshold.
       One row per lead_snp, ld_snp pair.

        * inputs:
            - clumped_file: plink output from --clump. contains one row per lead snps with snps in ld
                - full path
            - input_r2_file: found in the clump_snps.py output. contains lead and ld snp pair with r2 value
                - full path
            - upperinclusivebound and lowerexclusivebound
                - thresholds for r2. note lower bound is exclusive.

        * output:
            - input_df: see details below
            - one row for each input snp in ld with the lead snp

        * details:
            - locus_lead_snp:   lead gwas snp that defines the locus
            - input_snp:        gwas snps in ld with locus_lead_snp.
                                note:the locus_lead_snp is also present in the input_snp column with r2 == 1
            - lead_snp:         True/False on whether the input snp is a lead snp

            - note: if a lead snp has no ld parterns, then there is one row with lead snp with itself'''

    # create lead_snp, snp_in_ld_with_lead_snp pair (one per row)
    # created based on running plink --clump

    store_lead_snp = []
    store_input_snp = []

    # pdb.set_trace()


    clump_df = pd.read_csv(clumped_file, sep="\s+")

    for _, row in clump_df.iterrows():
        lead_snp = row.SNP
        sp2_row = row.SP2

        # each lead snp is also counted as a input_snp
        store_lead_snp.append(lead_snp)
        store_input_snp.append(lead_snp)

        if sp2_row == 'NONE':
            continue
        else:
            snps = [x[:-3] for x in sp2_row.split(",")]  # remove '(1)' at the end of each snp

            for ld_snp in snps:
                store_lead_snp.append(lead_snp)
                store_input_snp.append(ld_snp)

    # load gwas lead and ld snps
    rawinput_df = pd.DataFrame({'locus_lead_snp': store_lead_snp, 'input_snp': store_input_snp})
    r2_df = pd.read_csv(input_r2_file, sep="\t")
    r2_df.rename(columns={'ld_snp': 'input_snp', 'lead_snp': 'locus_lead_snp'}, inplace=True)

    # append r2
    merged_df = pd.merge(rawinput_df, r2_df, on=['locus_lead_snp', 'input_snp'], how='left')
    merged_df.loc[merged_df['locus_lead_snp'] == merged_df['input_snp'], 'R2'] = 1.0

    # filter on ld bin threshold
    # upperinclusivebound = 1.0
    # lowerexclusivebound = 0.9

    input_df = merged_df.loc[(merged_df['R2'] <= upperinclusivebound) & (merged_df['R2'] > lowerexclusivebound), ['locus_lead_snp', 'input_snp']].copy()

    input_df['lead_snp'] = input_df.loc[:, 'locus_lead_snp'] == input_df.loc[:, 'input_snp']
    input_df.reset_index(drop=True, inplace=True)
    return input_df


def combine_control_sets(control_sets_dir, clumped_dir, output_dir):

    tstart = time.time()

    logger.info("Combining all control sets into one file.")


    # SET UP FILES AND DIRS
    clumped_file = os.path.join(clumped_dir, "gwas_snps_clumped.tsv")
    input_r2_file = os.path.join(clumped_dir, "lead_gwas_snp_and_ld_snp_r2.tsv")
    matched_dir = os.path.join(control_sets_dir, 'matching_sets')

    # OUTPUT
    combined_sets_file_out = os.path.join(output_dir, 'combined_control_sets.tsv')



    ##
    ##   CREATE DF OF LEAD, INPUT SNP PAIR PER ROW
    ##

    # need this to create a table where each row is the lead snp or a snp in ld with the lead snp
    input_df = create_input_df(clumped_file, input_r2_file, upperinclusivebound=1.0, lowerexclusivebound=0.9)


    ##
    ##   PARSE CONTROL SET AND LINE UP WITH INPUT_DF
    ##

    # preallocate df with num columns equal to num control sets
    set_file_list = glob.glob(matched_dir+"/Set_*-ld_matched.txt")
    setnames = [os.path.split(x)[1].split("-")[0] for x in set_file_list]

    preallocate = pd.DataFrame(np.zeros([input_df.shape[0], len(set_file_list)])*np.nan, columns=setnames)
    input_df = pd.concat([input_df, preallocate], axis=1)


    ##
    ##   ADD EACH CONTROL SET TO A COLUMN
    ##

    store_mismatches = []
    for index, set_file in enumerate(set_file_list):

        setnum = os.path.split(set_file)[1]
        logger.debug("Lining up {}, out of {}".format(index, len(set_file_list)))
        input_df, mismatches = add_control_column(set_file, input_df)

        store_mismatches.append(mismatches)


    # check for any mismatches
    if len(np.hstack(store_mismatches)) > 0:
        logger.info("There were mismatches!")
        logger.info(np.hstack(store_mismatches))




    ##
    #   QC CONTROL SETS
    ##



    missing_controls_thresh = 0.05
    # number of missing snps / number of total input snps for each control set
    percent_missing = pd.isnull(input_df).sum()/input_df.shape[0]
    control_sets_to_remove = percent_missing[percent_missing > missing_controls_thresh].index
    clean_control_df = input_df.drop(control_sets_to_remove, axis=1)



    clean_control_df = input_df.copy()

    logger.info("Control sets were removed if they were missing more than {}% of required control snps.".format(missing_controls_thresh*100))
    logger.info("This resulted in removing {} sets, leaving {} control sets.".format(len(control_sets_to_remove), (clean_control_df.shape[1]-3)))

    ##
    #   WRITE COMBINED DATA
    ##

    clean_control_df.to_csv(combined_sets_file_out, sep="\t", index=False, na_rep="NA")
    logger.info("Filtered combined control set written to:{}".format(combined_sets_file_out, len(set_file_list)))
    logger.info("Done combining control sets into one file. It took {:.2f} minutes".format((time.time() - tstart)/60))

    return combined_sets_file_out


# -----------
# MAIN
# -----------

if __name__ == "__main__":


    control_sets_dir, clumped_dir = parse_input_args()
    combined_sets_file_out = combine_control_sets(control_sets_dir, clumped_dir)

