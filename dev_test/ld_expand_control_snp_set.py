#!/bin/python
# Take gwas clumped input snps and its corresponding SNPSNAP control snps and add control snps in ld with lead snps.
#       * Works on one control set per run of this script.
#
#
#       DEPENDENCIES:
#           * requires matched_snps.txt from SNPsnap output
#           * requires table with number of LD partners per LD bin for each lead snp
#
# Abin Abraham
# created on: 2018-08-31 00:51:40
# Python 3.6.1
# PLINK v1.90b3s

#######
# THIS SCRIPT WILL TRY TO EXPAND ALL CONTROL SETS



import os
import sys
import time
import glob
import argparse



import numpy as np
import pandas as pd



sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline")
from helper_get_ld_partners import get_ld_partners, get_ld_partners_modified, quantify_duplicates  # pylint: disable=E0401
from helper_general import Outputs, error_check_plink_run, warning_check_plink_run

import logging
logger = logging.getLogger('main.{}'.format(__name__))



# -----------
# PATHS
# -----------

THOU_GENOME_DIR = '1kg'
thous_gen_file = os.path.join(THOU_GENOME_DIR, "EUR.{}.phase3.nodups")



# -----------
# FUNCTIONS
# -----------

def parse_input_args():


    # TODO: delete (only for dev purposes)
    if len(sys.argv) != 11:

        matched_snps_file = "test_matched_snps/matched_snps.tsv"
        num_ld_snps_required_file = "test_clump/lead_gwas_snps_with_ldsnps_counts.txt"
        ld_for_control_snps_file = "ld_snps_for_control_snps"
        output_root = os.getcwd()
        contron_set_index = 1

    else:

        # TODO: convert all args to requried args
        parser = argparse.ArgumentParser(description='Run expand control set.')

        parser.add_argument('matched_snps_file',
                            action='store', type=str,
                            help="snpsnap matched_snps.txt file path")

        parser.add_argument('ldsnps_file',
                            action='store', type=str,
                            help='required ld snps per lead snp; output from clump_snps.py')

        parser.add_argument('ld_snps_for_control_snps_dir',
                            action='store', type=str,
                            help='directory to plink output from ld expanding control snps')

        parser.add_argument('set_index',
                            action='store', type=int,
                            help='set index for snpsnap matched snps, start counting at 1')

        parser.add_argument('analysis_name',
                            action='store', type=str,
                            help="the name of this analysis")


        parser.add_argument('output_dir', dest='output_dir',
                            default="control_ld_snps_hou",
                            action='store', type=str,
                            help='ouput_directory')

        parser.add_argument('analysis_name',
                            action='store', type=str,
                            help="the name of this analysis")

        # retrieve passed arguments
        args = parser.parse_args()
        matched_snps_file = args.matched_snps_file
        num_ld_snps_required_file = args.ldsnps_file
        ld_for_control_snps_file = args.ld_snps_for_control_snps_dir
        output_root = args.output_dir
        contron_set_index = args.set_index
        analysis_name = args.analysis_name

    return matched_snps_file, num_ld_snps_required_file, ld_for_control_snps_file, output_root, contron_set_index, analysis_name


def set_up_outputs(OutObj, set_index):

    # set up ouput dirs
    OutObj.add_output('matched_controlset_dir', "matching_sets", mkdir=True, add_root=True)
    OutObj.add_output('matching_summary_dir', "matching_summary", mkdir=True, add_root=True)


    # set up ouput files
    OutObj.add_output('ld_matched_output_file', "Set_{}-ld_matched.txt".format(set_index), custom_root=OutObj.get('matched_controlset_dir'))
    OutObj.add_output('matching_quality_file',"Set_{}-matching_quality.tsv".format(set_index), custom_root=OutObj.get('matching_summary_dir'))


    return OutObj


def ld_expand_control_snps(matched_snps_file, required_ldsnps_file, plink_ld_output_dir, output_root, analysis_name, set_index):
    tstart = time.time()

    logger.info("Starting to LD expand control snps for set {}.".format(set_index))
    output_dir = os.path.join(output_root, '{}_ld_expand_control_snps'.format(analysis_name))
    OutObj = Outputs(output_dir, overwrite=True)

    OutObj = set_up_outputs(OutObj, set_index)

    # load table with num ld control snps for each input snp
    required_df = pd.read_csv(required_ldsnps_file, sep="\t")


    # load snpsnap control snps
    snpsnap_df = pd.read_csv(matched_snps_file, sep="\t")
    set_label = snpsnap_df.columns[set_index]


    # TODO: harmonize final snp list to match output from snpsnap
    #

    # run plink ld and randomly pick ld controls
    mstart = time.time()
    count_dict = dict()
    matched_results_dict = dict()
    counter = 0


    for chrom in np.arange(1,23,1):

        # load one control set
        temp_chrm = snpsnap_df.loc[:, 'Set_{}'.format(set_index)].str.split(":", expand=True)
        chrom_bool = (temp_chrm.iloc[:,0]==str(chrom)) # 0 is indexing the chromosome column...

        # if there are no control snps on this chromosome.
        if not chrom_bool.any():
            continue

        # make a dictionary of input_snp:control_snp
        input_lead_snps,this_control_set = zip(*snpsnap_df.loc[chrom_bool, ['snps_to_match', 'Set_{}'.format(set_index)]].values)
        input_control_dict = dict(zip(input_lead_snps, this_control_set))



        # get all the SNPs in LD with the control SNP
        # PLINK_OUT_FILE contains all snps in ld with the control snps
        PLINK_OUT_FILE=os.path.join(plink_ld_output_dir,"ld_for_control_snps_chr{}")
        ld_df = pd.read_csv(PLINK_OUT_FILE.format(chrom)+".ld", sep="\s+")

        store_matched_snps = pd.DataFrame({'gwas_snp': [], 'control_snp': [], 'lead_snp': [], 'r2': []})
        for input_snp, control_snp in input_control_dict.items():

            logger.debug("Finding LD buddies for snp #{} of {}".format(counter+1, snpsnap_df.shape[0]))
            filt_ld_df = ld_df.loc[ld_df['SNP_A'] == control_snp].copy()

            # get ld controls
            # required_df holds the number of ld snps required for specified ld level in ld_bins_to_match
            # instructions for how to set ld_bins_to_match are in the function documentation
            matched_df, counts, flags = get_ld_partners_modified(filt_ld_df, required_df, input_snp, control_snp,  ld_bins_to_match=['ld<=1.0'])
            matched_results_dict[input_snp] = matched_df.to_dict()

            ld_score = matched_df.r2.sum()
            all_partners = matched_df.iloc[:, 1]
            # evaluate duplicate ld controls
            frac_non_redund = quantify_duplicates(all_partners)

            # unpack matching metrics
            count_dict[input_snp] ={"control_snp": [control_snp],
                                    "required_ld_snps": [counts['total_req_snps']],
                                    "ld_snps_matched": [counts['num_matched_snp']],
                                    "ld_score": [ld_score],
                                    "replacements": [counts['replacement_count']],
                                    "missing_ldbins": [counts['imperfect_match_count']],
                                    "no_ld_snps_available": [flags['no_ld_snps_available_flag']],
                                    "frac_non_redund": [frac_non_redund]}

            counter += 1

    logger.debug("ld expanding done: {:.4f} seconds".format( (time.time() - mstart)))

    # unwrap ld matched snps
    matched_df = pd.DataFrame()
    for index, inner_dict in matched_results_dict.items():
        matched_df = matched_df.append(pd.DataFrame(inner_dict))


    counts_df = pd.DataFrame.from_dict(count_dict, "index").applymap(lambda x: x[0])
    counts_df.reset_index(inplace=True)
    counts_df.rename(columns={'index':'gwas_input_snp'}, inplace=True)

    ###
    # create summary score
    ###

    # add the original control snp
    total_matched_snps = np.sum(counts_df.ld_snps_matched.apply(lambda x: x+1))  # add one to include original control snp
    counts_df['wt_frac_nr'] = (counts_df.frac_non_redund *
                               counts_df.ld_snps_matched.apply(lambda x: x+1))/total_matched_snps
    counts_df['frac_matched'] = (counts_df.ld_snps_matched/counts_df.required_ld_snps)
    counts_df.loc[counts_df['required_ld_snps'] == 0, 'frac_matched'] = 1


    final_score_frac_non_redundant = np.sum(counts_df.wt_frac_nr)
    final_score_mean_frac_matched = np.mean(counts_df.frac_matched)


    ###
    # write summary
    ###

    # write summary of matching quality
    new_col_order = ['gwas_input_snp', 'control_snp', 'required_ld_snps',
                     'ld_snps_matched', 'frac_matched', "ld_score", 'no_ld_snps_available',
                     'missing_ldbins', 'replacements', 'frac_non_redund',
                     'wt_frac_nr']



    counts_df = counts_df.loc[:, new_col_order].copy()
    counts_df.sort_values('frac_matched', inplace=True, ascending=False)
    counts_df.to_csv(OutObj.get('matching_quality_file'), sep="\t", index=False)


    # write/save ld control set
    matched_df.to_csv(OutObj.get('ld_matched_output_file'), sep="\t", index=False)

    logger.info("Done getting ld controls for set {}, it took {:.2f} minutes for {} input SNPs.".format(set_index,
        (time.time()-tstart)/60, snpsnap_df.shape[0]))

    logger.debug("*** MEAN MATCH FRACTION: {:.3f}".format(final_score_mean_frac_matched))
    logger.debug("*** WEIGHTED  FRACTION NON-REDUNDANT: {:.3f}".format(final_score_frac_non_redundant))

    return OutObj


# -----------
# MAIN
# -----------
if __name__ == "__main__":

    matched_snps_file, required_ldsnps_file,plink_ld_output_dir, output_root, set_index, analysis_name = parse_input_args()

    OutputObj = ld_expand_control_snps(matched_snps_file, required_ldsnps_file,plink_ld_output_dir, output_root, set_index, analysis_name)

