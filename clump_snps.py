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
import sys
import glob
import time
import argparse

import numpy as np
import pandas as pd
from subprocess import Popen, PIPE

from helper_calc_r2 import calc_r2_for_input_snps #pylint:E0401
from helper_general import Outputs, error_check_plink_run, warning_check_plink_run
from helper_clump_snps import write_gwas_snps_by_chr, get_list_of_ld_snps, get_r2_for_lead_ld_snps, bin_ldsnp_per_leadsnp

tstart = time.time()

from datetime import datetime
DATE = datetime.now().strftime('%Y-%m-%d')


import logging
logger = logging.getLogger('main.{}'.format(__name__))



# -----------
# DEPENDENCIES
# -----------

# TODO: maybe centralize all the required paths?
# THOU_GENOME_DIR = "/dors/capra_lab/projects/gwas_allele_age_evolution/data/1KG_snpsnap/nodups"
# THOU_GENOME_DIR = '1kg'
# thous_gen_file = os.path.join(THOU_GENOME_DIR, "EUR.{}.phase3.nodups")


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

        parser.add_argument('-n','--analysis_name', dest='analysis_name',
                            action='store', type=str,
                            help='the name of this analysis')

        argparser = parser.parse_args()
        gwas_summary_file = argparser.gwas_file
        thous_gen_file = argparser.one_kg_file
        output_root = argparser.output_root
        analysis_name = argparser.analysis_name


    return gwas_summary_file,thous_gen_file,  output_root, analysis_name

def set_up_outputs(OutputObj):

    # set up ouput dirs
    OutputObj.add_output('gstats_by_chr_dir', "gwas_stats_by_chr", mkdir=True, add_root=True)
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

def clump_snps(gwas_summary_file, output_root, analysis_name, thous_gen_file):

    output_dir = os.path.join(output_root, '{}_clump'.format(analysis_name))
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)

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
    write_gwas_snps_by_chr(gwas_df, OutObj.get('gstats_by_chr_dir'), OutObj.get('gstats_pos_by_chr_dir'), gsum_header['chr'], gsum_header['basepair'], gwasfilename, gsum_header['pvalue'])


    ##
    ## run plink clump by chromosome
    ##

    logger.info(f"Running plink clump on gwas variants by chromosome.")
    keep_warnings = []
    for thisfile in glob.glob(OutObj.get('gstats_by_chr_dir')+"/chr*"):

        # inputs and outputs for plink call
        this_chr = os.path.split(thisfile)[1].split("_")[0]  # e.g. "chr3"
        chr_num = this_chr.split('chr')[1]
        outputfile = os.path.join(OutObj.get('plink_clump_output_dir'), this_chr)


        plinkcmd = f'plink --bfile {thous_gen_file.format(chr_num)} --clump {thisfile} --clump-p1 0.000001 --clump-p2 0.000001 --out {outputfile}'

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


# -----------
# MAIN
# -----------

if __name__ == "__main__":


    gwas_summary_file, thous_gen_file,  output_dir, analysis_name = parse_input_args()


    _ = clump_snps(gwas_summary_file,  output_dir, analysis_name, thous_gen_file)
