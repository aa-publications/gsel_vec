#!/bin/python
# Take GWAS summary statistics file -->
#               1) run plink clump
#               2) bin all gwas variants by different levels of LD
#
#
# DEPENDENCIES
#       1) check whether the p-value in gwas input file is in scientific e notation only!
#
#
#
#
# Abin Abraham
# created on: 2018-08-26 15:14:16
# modified on: 2020-01-05 08:16:24
# plink 1.90b
# python 3.6


import os
import re
import sys
import glob
import time
import argparse

import numpy as np
import pandas as pd
from subprocess import Popen, PIPE


from helper_calc_r2 import calc_r2_for_input_snps, ld_expand_snp_list #pylint:E0401
from helper_general import Outputs, error_check_plink_run, warning_check_plink_run
from helper_clump_snps import write_gwas_snps_by_chr, get_list_of_ld_snps, get_r2_for_lead_ld_snps, bin_ldsnp_per_leadsnp, write_snp_list_by_chr

tstart = time.time()

from datetime import datetime
DATE = datetime.now().strftime('%Y-%m-%d')


import logging
logger = logging.getLogger('main.{}'.format(__name__))

# %%


# -----------
# DEPENDENCIES
# -----------


# -----------
# FUNCION
# -----------


def parse_log_for_err(plink_clump_output_dir, missing_snps_file, plink_log_errors_file):

    store_variants = list()
    store_errors = list()

    for thislog in glob.glob(plink_clump_output_dir+"/chr*.log"):
        with open(thislog, 'r') as fh_log:

            for line in fh_log:
                if line.startswith('Warning:') and line.endswith('top variant.\n'):
                    missing_variant = line.split()[1][1:-1]
                    store_variants.append(missing_variant)

                if line.startswith("Error:"):
                    store_errors.append("{}: {}".format(thislog, line))

    if store_variants:
        with open(missing_snps_file, 'w') as fout:
            fout.write("\n".join(store_variants))

    if store_errors:
        with open(plink_log_errors_file, 'w') as ferr:
            ferr.write("\n".join(store_errors))

    return store_errors, store_variants

def concat_plink_clump_output(plink_clump_output_dir):
    store_clump_df = pd.DataFrame()
    for thisclump in glob.glob(plink_clump_output_dir+"/chr*.clumped"):

        clump_df = pd.read_csv(thisclump, sep="\s+")
        store_clump_df = store_clump_df.append(clump_df)

    # format and write concatenated clump file
    store_clump_df['CHR_BP'] = store_clump_df['CHR'].map(str) + ":" + store_clump_df['BP'].map(str)

    return store_clump_df

def parse_input_args():


    # TODO: delete (only for dev purposes)
    if len(sys.argv) != 5:
        print("running dev args...delete later...")
        gwas_summary_file = "test_data/input_gwas_file.tsv"
        output_root= os.getcwd()
    else:


        parser = argparse.ArgumentParser(description='Will LD clump from GWAS summary stats.')

        parser.add_argument('-g', '--gwas_summary_stat_file', dest='gwas_file',
                            action='store', type=str,
                            help="path to gwas summary file")

        parser.add_argument('-k', '--1kg_file', dest='one_kg_file',
                            action='store', type=str,
                            help="path to gwas summary file")

        parser.add_argument('-o','--output', dest='output_root',
                            action='store', type=str,
                            help='path to dir where output will be stored')


        argparser = parser.parse_args()
        gwas_summary_file = argparser.gwas_file
        thous_gen_file = argparser.one_kg_file
        output_root = argparser.output_root



    return gwas_summary_file,thous_gen_file,  output_root

def set_up_outputs(OutputObj):

    # set up ouput dirs
    OutputObj.add_output('by_chr_dir', "gwas_stats_by_chr", mkdir=True, add_root=True)
    OutputObj.add_output('gstats_pos_by_chr_dir', "gwas_stats_pos_only_by_chr", mkdir=True, add_root=True)
    OutputObj.add_output('plink_clump_output_dir', "plink_clump_output", mkdir=True, add_root=True)
    OutputObj.add_output('gstats_r2_dir', "pairwise_r2_between_gwas_snps", mkdir=True, add_root=True)

    # set up ouput files
    OutputObj.add_output('reprint_gwas_file', "input_gwas_file.tsv", add_root=True)
    OutputObj.add_output('missing_snps_file', "input_gwas_snps_not_found.txt", add_root=True)
    OutputObj.add_output('lead_snps_file', "lead_gwas_snps.tsv", add_root=True)
    OutputObj.add_output('plink_log_errors_file', "plink_clump_errors.txt", add_root=True)
    OutputObj.add_output('clumped_snps_file', "gwas_snps_clumped.tsv", add_root=True)
    OutputObj.add_output('bin_by_ld_file', "lead_gwas_snps_with_ldsnps_counts.txt", add_root=True)
    OutputObj.add_output('lead_snp_ld_pairs_r2', "lead_gwas_snp_and_ld_snp_r2.tsv", add_root=True)
    OutputObj.add_output('log_file', 'log_clump_snps_{}.txt'.format(DATE),add_root=True)


    return OutputObj


def set_up_outputs_for_snp_list(OutObj):
    OutObj.add_output('by_chr_dir', "snps_by_chr", mkdir=True, add_root=True)

    # OutObj.add_output('ld_expand_r2_dir', "ld_expand_snp_list", mkdir=True, add_root=True)
    OutObj.add_output('plink_r2_output', "ld_expand_plink_output", mkdir=True, add_root=True)
    OutObj.add_output('all_input_snps_ld_expanded', "all_input_snps_ld_expanded.tsv", add_root=True)


    OutObj.add_output('pairwise_r2_of_input_snps', "pairwise_r2_of_input_snps", add_root=True, mkdir=True)
    OutObj.add_output('bin_by_ld_file', "input_snps_with_binned_ldsnps_counts.tsv", add_root=True)
    OutObj.add_output('lead_snp_ld_pairs_r2', "input_and_ld_snps.tsv", add_root=True)
    OutObj.add_output('non_indep_snps_file', "non_indep_snps.tsv", add_root=True)

    # set up ouput files
    OutObj.add_output('missing_snps_file', "input_snps_not_in_1kg_excluded.txt", add_root=True)
    OutObj.add_output('plink_log_errors_file', "plink_ld_expand_errors.txt", add_root=True)

    return OutObj

def force_input_snp_in_first_col(keep_autosomal_snps, og_store_ld_df):
    """ ensure that the input snp is in SNP_A column; if there are no LD snps, then there is a row for that SNP with 'NONE' as placeholders """

    store_ld_df = og_store_ld_df.copy()

    # check that all input snps are in SNP_A column
    remaining_autosomal_snps_in_SNP_A = set(keep_autosomal_snps).difference(set(store_ld_df.SNP_A.unique()))
    input_autosomal_snps_in_SNP_B = set(keep_autosomal_snps).difference(set(store_ld_df.SNP_B.unique())).difference(remaining_autosomal_snps_in_SNP_A) # remove any snps found in SNP_A column already


    if len(remaining_autosomal_snps_in_SNP_A) != 0:

        # chceck if any missing input snps from SNP_A are in SNP_B column
        remaining_input_snps_in_snp_b = set(remaining_autosomal_snps_in_SNP_A).difference(set(store_ld_df.SNP_B.unique()))

        if len(remaining_input_snps_in_snp_b) > 0:

            # SWAP THE SNP_A and SNP_B labels
            temp_snpb = store_ld_df.loc[store_ld_df['SNP_B'].isin(remaining_input_snps_in_snp_b)].copy()
            assert (temp_snpb.columns.values.tolist() == ['SNP_A', 'SNP_B', 'R2', 'snpA_B', 'snpB_A']), 'columns are not what is expected'
            temp_snpb.columns = ['SNP_B','SNP_A','R2','snpB_A','snpA_B']

            # drop these snps from the initial df and add the swapped columsn back in
            store_ld_df.drop(temp_snpb.index, axis=0, inplace=True)

            # add the swapped snps back in
            store_ld_df = store_ld_df.append(temp_snpb, sort=True)

        # check again if there are missing snps from SNP_A column
        remaining_autosomal_snps_in_SNP_A_after_swap = set(keep_autosomal_snps).difference(set(store_ld_df.SNP_A.unique()))



        if len(remaining_autosomal_snps_in_SNP_A_after_swap) > 0:


            add_miss_snp_df =  pd.DataFrame({'SNP_A':list(remaining_autosomal_snps_in_SNP_A_after_swap),'SNP_B':"NONE",'R2':"NONE",'snpA_B':"NONE",'snpB_A':"NONE"})
            store_ld_df = store_ld_df.append(add_miss_snp_df, sort=True)


    # reorder columns
    store_ld_df = store_ld_df.loc[:, ['SNP_A','SNP_B','snpA_B','snpB_A','R2']]
    return store_ld_df

def clump_snps(gwas_summary_file, output_root, thous_gen_file, lead_snp_min_gwas_pvalue=0.00000005, ld_snps_min_gwas_pvalue=0.00000005, min_r2_to_clump=0.9, min_kb_from_index_to_clump=250):
    
    

    
    output_dir = os.path.join(output_root, 'ld_clump_inputs')
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)


    # TODO: figure out how to deal with input data...
    # excpected header names for gwas summary statistics file
    gsum_header = {'rsID':"snp",
                    'chr':"chr",
                    'basepair':"pos",
                    'pvalue':"p"}



    ##
    ## parse gwas input snps
    ##

    gwas_df = pd.read_csv(gwas_summary_file, sep="\t")
    gwas_df.to_csv(OutObj.get('reprint_gwas_file'), sep="\t", index=False)  # write a copy of the input data
    logger.info(f"Loaded GWAS summary stats with {gwas_df.shape[0]} rows.")

    # write input snps by chromosome
    logger.info(f"Splitting GWAS summary stats by chromosome.")
    gwasfilename = os.path.split(gwas_summary_file)[1]
    write_gwas_snps_by_chr(gwas_df, OutObj.get('by_chr_dir'), OutObj.get('gstats_pos_by_chr_dir'), gsum_header['chr'], gsum_header['basepair'], gwasfilename, gsum_header['pvalue'])


    ##
    ## run plink clump by chromosome
    ##

    logger.info(f"Running plink clump on gwas variants by chromosome.")
    keep_warnings = []
    for thisfile in glob.glob(OutObj.get('by_chr_dir')+"/chr*"):

        # inputs and outputs for plink call
        this_chr = os.path.split(thisfile)[1].split("_")[0]  # e.g. "chr3"
        chr_num = this_chr.split('chr')[1]
        outputfile = os.path.join(OutObj.get('plink_clump_output_dir'), this_chr)

    
        plinkcmd = f'plink --bfile {thous_gen_file.format(chr_num)} --clump {thisfile} --clump-kb {min_kb_from_index_to_clump} --clump-r2 {min_r2_to_clump} --clump-p1 {lead_snp_min_gwas_pvalue} --clump-p2 {ld_snps_min_gwas_pvalue} --out {outputfile}'
        logger.debug(plinkcmd)

        # run plink cmd
        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()

        _ = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        keep_warnings.append(any_warnings)

    if np.any(keep_warnings):
        print("Plink threw warnings...check the log file.")


    if (len(glob.glob(OutObj.get('plink_clump_output_dir')+"/chr*.clumped")) == 0):
        raise RuntimeError(f'No plink clumped results in:\n{OutObj.get("plink_clump_output_dir")}')
    else:
        _,store_variants = parse_log_for_err(OutObj.get('plink_clump_output_dir'), OutObj.get('missing_snps_file'), OutObj.get('plink_log_errors_file'))

    # write all plink clump output by chromosome to one file
    store_clump_df = concat_plink_clump_output(OutObj.get('plink_clump_output_dir'))
    store_clump_df.rename(columns={'SNP':'lead_snp'}, inplace=True)
    store_clump_df.to_csv(OutObj.get('clumped_snps_file'), index=False, header=True, sep="\t")
    logger.debug("Combined all plink clump files by chromosome to: {}".format((OutObj.get('clumped_snps_file'))))


    ##
    ## calc r2 between gwas snps
    ##

    # make data frame with lead and ld snp with its r2
    logger.info("Calculating pairwise r2 for gwas variants.")

    r2_df = calc_r2_for_input_snps(OutObj.get('gstats_pos_by_chr_dir'), OutObj.get('gstats_r2_dir'), thous_gen_file)


    # load all lead_snp and ld_snp pairs
    lead_ld_df = pd.DataFrame()
    for this_chr_clump in glob.glob(OutObj.get("plink_clump_output_dir")+"/chr*.clumped"):

        temp_df  = pd.DataFrame(get_list_of_ld_snps(this_chr_clump), columns=['lead_snp','ld_snp'])
        lead_ld_df = lead_ld_df.append(temp_df)



    # get r2 for all lead and ld snp pairs
    store_lead_snp_ld_df = get_r2_for_lead_ld_snps(lead_ld_df, r2_df)
    store_ld_bins_df = bin_ldsnp_per_leadsnp(store_lead_snp_ld_df)
    store_ld_bins_df.reset_index(drop=False, inplace=True)
    store_ld_bins_df.rename(columns={'index':'lead_snp'}, inplace=True)




    ##
    ## write outputs
    ##

    # lead snps only file
    # TODO: delete the label to this file from the OUTOBJ object
    # lead_ld_df[~lead_ld_df.duplicated('lead_snp', keep='first')].to_csv(OutObj.get('lead_snps_file'), sep="\t", index=False, header=False, columns=['lead_snp'])
    # logger.info(f"Wrote lead snps only: {OutObj.get('lead_snps_file')}")

    # ld_bins
    store_ld_bins_df.to_csv(OutObj.get('bin_by_ld_file'), sep="\t", index=False)
    logger.debug(f"Wrote table of number of LD snps (binned): {OutObj.get('bin_by_ld_file')}")

    # lead snp, LD snp, r2 table
    store_lead_snp_ld_df.to_csv(OutObj.get('lead_snp_ld_pairs_r2'), sep="\t", index=False)
    logger.debug(f"Wrote table of lead and ld snp pairs with r2: {OutObj.get('lead_snp_ld_pairs_r2')}")

    logger.info(f"* While running plink clump, {len(store_variants)} input SNPs not found in 1KG written to:\n\t{OutObj.get('missing_snps_file')}.")
    logger.info(f"Done clumping GWAS summary stats. Found {store_ld_bins_df.lead_snp.nunique() } lead snps and took {(time.time()-tstart)/60:.2f} minutes.")

    return OutObj


def clump_snp_list(snps_list_file, output_root, thous_gen_file, min_r2_to_clump=0.9):


    # set up outputs
    output_dir = os.path.join(output_root, 'ld_clump_inputs')
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs_for_snp_list(OutObj)



    ###
    ###    load list of snp and check input format
    ###

    # laod input snps
    with open(snps_list_file, 'r') as fo:
        snp_file_contents = fo.read().splitlines()

    logger.info(f"Loaded {len(snp_file_contents):,} snps for analysis.")
    # ensure input snps follow the following format chr:pos
    assert np.all([ True if re.match(r"[1-9][0-9]*:[0-9]+", str(snp)) else False for snp in snp_file_contents ]), "check input snps are in formatted correctly (2:22222, chr:pos, without any prefix for chromosome)"

    # keep only autosomal snps
    keep_autosomal_snps = [snp  for snp in snp_file_contents if (int(snp.split(":")[0]) < 23)]
    non_auto_snps_removed = set(snp_file_contents).difference(set(keep_autosomal_snps))
    logger.info(f"Removed {len(non_auto_snps_removed)} snps for not being autosomal snps.")


    if len(keep_autosomal_snps) == 0:
        sys.exit("Error, script terminated because no autosomal snps left in the input snp list.")




    ###
    ###    ld expand each input snp
    ###

    # split by chromosome
    snp_list_by_chr_files_written = write_snp_list_by_chr(keep_autosomal_snps, OutObj.get('by_chr_dir'))


    # ld expand
    plink_r2_output_dir = OutObj.get('plink_r2_output')
    all_snps_ld_expanded_file = OutObj.get('all_input_snps_ld_expanded')
    store_ld_df = ld_expand_snp_list(snp_list_by_chr_files_written, plink_r2_output_dir, thous_gen_file)

    _,store_miss_variants  = parse_log_for_err(plink_r2_output_dir, OutObj.get('missing_snps_file'), OutObj.get('plink_log_errors_file'))



    ###
    ###    bin snps by ld
    ###

    # remove snps that were not found in 1kg database, so that these SNPs don't get filled in as having no SNP in LD with it
    keep_autosomal_snps = set(keep_autosomal_snps).difference(set(store_miss_variants))


    # force input snp in SNP_A and bin
    for_ld_binning_df = force_input_snp_in_first_col(keep_autosomal_snps, store_ld_df)
    for_ld_binning_df.rename(columns={'SNP_A':'lead_snp'}, inplace=True)
    for_ld_binning_df.R2 = for_ld_binning_df.R2.round(2)
    store_ld_bins_df  = bin_ldsnp_per_leadsnp(for_ld_binning_df)
    store_ld_bins_df.reset_index(drop=False, inplace=True)
    store_ld_bins_df.rename(columns={'index':'lead_snp'}, inplace=True)


    ###
    ###    calculate ld between input snps and check for ld independence
    ###


    pairwise_r2_df = calc_r2_for_input_snps(OutObj.get('by_chr_dir'), OutObj.get('pairwise_r2_of_input_snps'), thous_gen_file, min_r2_threshold=0.7)

    # write a list of input snps are within the given r2 threshold
    forced_pairwise_r2_df = force_input_snp_in_first_col(keep_autosomal_snps, pairwise_r2_df)
    non_indep_input_snps_df = pairwise_r2_df.loc[pairwise_r2_df['R2']> min_r2_to_clump].copy()



    ###
    ###    write
    ###

    # write input snps and ld snps
    lead_ld_snp_df = for_ld_binning_df.loc[for_ld_binning_df['R2'] > min_r2_to_clump, ['lead_snp','SNP_B','R2']].copy()
    lead_ld_snp_df.rename(columns={'SNP_B':'ld_snp'},inplace=True)
    lead_ld_snp_df.R2 = lead_ld_snp_df.R2.round(2)
    lead_ld_snp_df.to_csv(OutObj.get('lead_snp_ld_pairs_r2'), sep="\t", index=False, header=True)
    logger.debug(f"Wrote input and ld snp to: {OutObj.get('lead_snp_ld_pairs_r2')}")


    # write snps binned ld
    store_ld_bins_df.to_csv(OutObj.get('bin_by_ld_file'), sep="\t", index=False)
    logger.debug(f"Wrote table of number of LD snps (binned): {OutObj.get('bin_by_ld_file')}")


    # write non_indep snps
    write_non_indep_df = non_indep_input_snps_df.loc[non_indep_input_snps_df['SNP_A'].isin(keep_autosomal_snps)].copy()

    if write_non_indep_df.shape[0] > 0:
        write_non_indep_df.to_csv(OutObj.get('non_indep_snps_file') , sep="\t", index=False, header=True)
        logger.debug(f"Wrote input snps that are non-indepdent at r2 {min_r2_to_clump} to: {OutObj.get('non_indep_snps_file')}")


    # write input snps excluded from analysis
    if store_miss_variants:
        with open(OutObj.get('missing_snps_file'), 'w') as fo:
            fo.write('\n'.join(store_miss_variants))

    return OutObj






# -----------
# MAIN
# -----------

if __name__ == "__main__":


    gwas_summary_file, thous_gen_file,  output_dir = parse_input_args()


    _ = clump_snps(gwas_summary_file,  output_dir, thous_gen_file)
