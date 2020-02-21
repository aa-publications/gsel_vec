#!/bin/python
# This script will ...
#
#
#       DEPENDENCIES
#               * assumes that the headers for file with number of LD snps requires has headers in this format: 'ld<=0.0', 'ld<=0.1'
# Abin Abraham
# created on: 2018-09-03 00:53:37

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE


tstart = time.time()


import logging
logger = logging.getLogger('main.{}'.format(__name__))



# -----------
# FUNCTIONS
# -----------

def get_r2mask(right_bound, left_bound, ld_df, num_snps_req):

    # this ld bin had no matches
    imperfect_match_count = 0

    r2_mask = (ld_df['R2'] <= right_bound) & (ld_df['R2'] > left_bound)

    # if there are no snps in ld for these bounds
    while np.any(r2_mask) == False:

        imperfect_match_count = 1

        # decrease the ld bin by 0.1 and try again
        right_bound -= 0.1
        left_bound -= 0.1

        # choose from all available snps if left bound < 0.5
        if left_bound < 0.5:
            r2_mask = r2_mask | True
        else:
            r2_mask = (ld_df['R2'] <= right_bound) & (ld_df['R2'] > left_bound)

    return r2_mask, imperfect_match_count

def get_ld_partners(snp_ld_file, ld_partners_binned_df, gwas_input_snp, control_snp, ld_bins_to_match=['ld<=0.0', 'ld<=0.1', 'ld<=0.2',
                                                                                                       'ld<=0.3', 'ld<=0.4', 'ld<=0.5', 'ld<=0.6',
                                                                                                       'ld<=0.7', 'ld<=0.8', 'ld<=0.9', 'ld<=1.0']):
    """ Randomly samples snps (ld partners) in ld with given control snp.

    Dependencies: assumes that the headers for ld_partners_binned_df follow this format: 'ld<=0.0', 'ld<=0.1', etc.

    If required ld partners is 0, return the control snp(CHROM:BP).
    Otherwies:
        1) if ld file exists randomly sample ld snps for a given ld bin
        2) if ld file is missing, return control snp.

    Corner Cases:
        N = number of snps required for LD bin.
        1) if there are no SNPs in LD with control SNP, return the control snp(CHROM:BP).
        2) if there are no SNPs for a particular LD BIN, then keep checking the next higher bin until a match is found.
            2a) if no matches are found after reach r2 == 1, then randomly sample from whatever SNPs are avabilable.
        3) if N > available SNPs to draw from for a particular LD bin, sample with replacement. """

    # set up
    # grab count of snps required for the specified LD bins
    # req_ld_count_array = ld_partners_binned_df.loc[ld_partners_binned_df['lead_snp'] == gwas_input_snp, ld_bins_to_match ].values.flatten()
    # ld_bins_to_match = ['ld<=0.9','ld<=1.0']
    req_ld_count_array = ld_partners_binned_df.loc[ld_partners_binned_df['lead_snp']
                                                   == gwas_input_snp, ld_bins_to_match]

    pos_ld_bins_index = np.where(req_ld_count_array.values.flatten())[0]
    pos_ld_bins_colnames = req_ld_count_array.columns[pos_ld_bins_index]
    total_req_snps = np.sum(req_ld_count_array.values.flatten())

    replacement_counter = 0
    missing_ld_file_flag = False
    no_ld_snps_req_flag = False
    no_ld_snps_available_flag = False
    imperfect_match_count = 0

    store_partners = [control_snp]
    store_r2 = [1]
    logger.debug("Getting LD partners for control snp {}".format(control_snp))
    # if no ld partners required
    if pos_ld_bins_index.size == 0:
        # logger.info("For control snp {}, no LD partners were required.\n".format(control_snp))
        no_ld_snps_req_flag = True

    # if ld file exists
    elif os.path.isfile(snp_ld_file):

        ld_df = pd.read_csv(snp_ld_file, sep="\s+")  # pylint: disable=W1401

        # remove self:self pair
        ld_df.drop(ld_df.loc[ld_df['BP_B'] == ld_df["BP_A"]].index, inplace=True)
        ld_df["CHR_B:BP_B"] = ld_df["CHR_B"].map(str) + ":" + ld_df["BP_B"].map(str)

        # no snps in ld with control snp
        if ld_df.shape[0] == 0:
            logger.warning("no snps in ld with control snpt to choose from...")
            no_ld_snps_available_flag = True
        else:

            # loop through each LD BIN
            for ld_bin in pos_ld_bins_colnames:

                # derive ld bin from ld_header
                right_inclusive_bound = float(ld_bin.split("=")[1])
                left_exclusive_bound = right_inclusive_bound - 0.1
                num_snps_req = req_ld_count_array[ld_bin].values[0]

                r2_masks, imperfect_match = get_r2mask(right_inclusive_bound, left_exclusive_bound, ld_df, num_snps_req)
                imperfect_match_count += imperfect_match

                # logger.debug("number of potential LD snps: {}".format(len(np.where(r2_masks)[0])))
                # logger.debug("number of required LD snps: {}".format(num_snps_req))

                if len(np.where(r2_masks == True)[0]) < num_snps_req:

                    replacement_counter += 1
                    ld_control_snps = ld_df.loc[r2_masks, ["CHR_B:BP_B", 'R2']].sample(
                        n=num_snps_req, replace=True, random_state=1).values

                else:
                    ld_control_snps = ld_df.loc[r2_masks, ["CHR_B:BP_B", 'R2']].sample(
                        n=num_snps_req, replace=False, random_state=1).values

                store_partners.append(ld_control_snps[:, 0].tolist())
                store_r2.append(ld_control_snps[:, 1].tolist())

    # if no ld file exists for control snp
    else:
        missing_ld_file_flag = True
        logger.info("No LD file for {} control snp was found!".format(control_snp))

    all_partners_array = np.hstack(store_partners)
    all_r2_array = np.hstack(store_r2)

    num_matched_snp = len(all_partners_array) - 1  # remove original control snp

    flags = {'missing_ld_file_flag': missing_ld_file_flag,
             'no_ld_snps_req_flag': no_ld_snps_req_flag,
             'no_ld_snps_available_flag': no_ld_snps_available_flag}

    counts = {'num_matched_snp': num_matched_snp,
              'total_req_snps': total_req_snps,
              'replacement_count': replacement_counter,
              'imperfect_match_count': imperfect_match_count}

    num_rows = len(all_partners_array)
    matched_df = pd.DataFrame({'gwas_snp': [gwas_input_snp]*num_rows,
                               'control_snp': all_partners_array,
                               'r2': all_r2_array,
                               'lead_snp': np.hstack([True, [False]*(num_rows-1)])})
    matched_df = matched_df.loc[:, ['gwas_snp', 'control_snp', 'lead_snp', 'r2']]

    # print summary of matching
    if missing_ld_file_flag:
        logger.warning("For control snp {}, no ld file. (Likely, snps is missing in plink data).\n".format(control_snp))
    elif no_ld_snps_available_flag:
        logger.warning("For control snp {}, there were no ld snps when we needed to draw ld snps.\n".format(control_snp))
    elif imperfect_match_count != 0:
        logger.info("For control snp {}, {} ld bins did not have matching snps.".format(
            control_snp, imperfect_match_count))
        logger.info("Found {} out of {} required ld controls, with {} replacements\n".format(
            num_matched_snp, total_req_snps, replacement_counter))
    else:
        logger.info("For control snp {}, found {} out of {} required ld controls.".format(
            control_snp, num_matched_snp, total_req_snps))
        logger.info("Sampled with {} replacement(s).\n".format(replacement_counter))

    return matched_df, counts,  flags

def get_ld_partners_modified(ld_df, ld_partners_binned_df, gwas_input_snp, control_snp, ld_bins_to_match=['ld<=0.0', 'ld<=0.1', 'ld<=0.2',
                                                                                                       'ld<=0.3', 'ld<=0.4', 'ld<=0.5', 'ld<=0.6',
                                                                                                       'ld<=0.7', 'ld<=0.8', 'ld<=0.9', 'ld<=1.0']):
    """ Randomly samples snps (ld partners) in ld with given control snp.

    IMPORTANT:
        * assumes that the headers for ld_partners_binned_df follow this format: 'ld<=0.0', 'ld<=0.1', etc.
        * if you want to match for LD where
            * 0.9 < r2 <= 1.0 ==> then ld_bins_to_match should only be 'ld<=1.0'
            * 0.8 < r2 <= 1.0 ==> then ld_bins_to_match should only be ['ld<=0.9', 'ld<=1.0']


    OTHER:

    If required ld partners is 0, return the control snp(CHROM:BP).
    Otherwise:
        1) if ld file exists randomly sample ld snps for a given ld bin
        2) if ld file is missing, return control snp.

    Corner Cases:
        N = number of snps required for LD bin.
        1) if there are no SNPs in LD with control SNP, return the control snp(CHROM:BP).
        2) if there are no SNPs for a particular LD BIN, then keep checking the next higher bin until a match is found.
            2a) if no matches are found after reach r2 == 1, then randomly sample from whatever SNPs are avabilable.
        3) if N > available SNPs to draw from for a particular LD bin, sample with replacement. """


    # pdb.set_trace()

    # get number of required snps
    req_ld_count_array = ld_partners_binned_df.loc[ld_partners_binned_df['lead_snp'] == gwas_input_snp, ld_bins_to_match]
    pos_ld_bins_index = np.where(req_ld_count_array.values.flatten())[0]
    pos_ld_bins_colnames = req_ld_count_array.columns[pos_ld_bins_index]
    total_req_snps = np.sum(req_ld_count_array.values.flatten())

    # initialize
    replacement_counter = 0
    no_ld_snps_req_flag = False
    no_ld_snps_available_flag = False
    imperfect_match_count = 0
    store_partners = [control_snp]
    store_r2 = [1]


    logger.debug("Getting LD partners for control snp {}".format(control_snp))


    if (total_req_snps == 0):

        no_ld_snps_req_flag = True


    else:

        # remove self:self pair
        ld_df.drop(ld_df.loc[ld_df['BP_B'] == ld_df["BP_A"]].index, inplace=True)
        ld_df["CHR_B:BP_B"] = ld_df["CHR_B"].map(str) + ":" + ld_df["BP_B"].map(str)


        if (ld_df.shape[0] == 0):

            # no snps in ld with control snp
            logger.debug("no snps in ld with control snp ({}) to choose from...".format(control_snp))
            no_ld_snps_available_flag = True


        else:

            # loop through each LD BIN
            for ld_bin in pos_ld_bins_colnames:

                # derive ld bin from ld_header
                right_inclusive_bound = float(ld_bin.split("=")[1])
                left_exclusive_bound = right_inclusive_bound - 0.1
                num_snps_req = req_ld_count_array[ld_bin].values[0]

                r2_masks, imperfect_match = get_r2mask(right_inclusive_bound, left_exclusive_bound, ld_df, num_snps_req)
                imperfect_match_count += imperfect_match

                logger.debug("number of potential LD snps: {}".format(len(np.where(r2_masks)[0])))
                logger.debug("number of required LD snps: {}".format(num_snps_req))

                if (len(np.where(r2_masks == True)[0]) < num_snps_req):

                    replacement_counter += 1
                    ld_control_snps = ld_df.loc[r2_masks, ["CHR_B:BP_B", 'R2']].sample(n=num_snps_req, replace=True, random_state=1).values

                else:
                    ld_control_snps = ld_df.loc[r2_masks, ["CHR_B:BP_B", 'R2']].sample(n=num_snps_req, replace=False, random_state=1).values

                store_partners.append(ld_control_snps[:, 0].tolist())
                store_r2.append(ld_control_snps[:, 1].tolist())


    all_partners_array = np.hstack(store_partners)
    all_r2_array = np.hstack(store_r2)
    num_matched_snp = len(all_partners_array) - 1  # remove original control snp

    flags = {'no_ld_snps_req_flag': no_ld_snps_req_flag,
             'no_ld_snps_available_flag': no_ld_snps_available_flag}

    counts = {'num_matched_snp': num_matched_snp,
              'total_req_snps': total_req_snps,
              'replacement_count': replacement_counter,
              'imperfect_match_count': imperfect_match_count}

    num_rows = len(all_partners_array)
    matched_df = pd.DataFrame({'gwas_snp': [gwas_input_snp]*num_rows,
                               'control_snp': all_partners_array,
                               'r2': all_r2_array,
                               'lead_snp': np.hstack([True, [False]*(num_rows-1)])})
    matched_df = matched_df.loc[:, ['gwas_snp', 'control_snp', 'lead_snp', 'r2']]


    return matched_df, counts,  flags


def quantify_duplicates(snps_array):
    ''' fraction of non-redundant values (:= 1 - fraction of all non-first-instance duplicates) '''

    _, uq_counts = np.unique(snps_array, return_counts=True)

    frac_redundant = np.sum(uq_counts - 1)/np.sum(uq_counts)
    frac_non_redund = 1 - frac_redundant

    return frac_non_redund
