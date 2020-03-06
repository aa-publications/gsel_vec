#!/bin/python
# Calculates pairwise ld for input gwas snps using plink.
# This is a helper function to be loaded with clump_snps.py.
#
#
#           NOTES:
#           All required files and directories must be passed in.
#           See function doc string for details.
#
#           DEPENDENCIES:
#           * requires PLINK v1.90b3s
#
#           OUTPUTS:
#           * will be written to directory; see doc string for details.
#
#
# Abin Abraham
# created on: 2018-08-26 15:14:16

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE


from helper_general import Outputs, safe_mkdir, error_check_plink_run, warning_check_plink_run


import logging

logger = logging.getLogger('main.{}'.format(__name__))


# -----------
# FUNCTIONS
# -----------



def calc_r2_for_input_snps(snp_list_by_chr_files_dict, gwas_snps_by_chr_files_dict, output_dir, thous_gen_file, min_r2_threshold=0.7):
    """"Calculates pariwise LD (r2) using plink --r2 inter-chr. See below for parameter details

        INPUTS:
            * snp_list_by_chr_files_dict: dictionary with [chromsome_number: path_to_file ] for files containing clumped gwas stats by chromosomes
            * gwas_snps_by_chr_files_dict:  dictionary with [chromsome_number: path_to_file ] for files containign all gwas variants by chromosome
            * output_dir: base directory where all resutls will be store; if it doesn't exist, it will create.
                > subdirectories will not be overwritten, will exit script if it already exists.

            * thous_gen_file: path to Thousand Genomes PLINK file; no further filtering is done


        OUTPUTS:
            * plink output is written to: plink_ld_output/.
            * pairwise ld is combined across chromosomes and written to: all_pairwise_ld.tsv
            * plink screen output is written to: plink_ld_stdout.txt"""

    ##
    ## setup
    ##

    tstart = time.time()

    # set up output folders and files
    ld_summary_file = os.path.join(output_dir, "all_pairwise_ld.tsv")
    plink_output_dir = os.path.join(output_dir, 'plink_ld_output')


    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    safe_mkdir(plink_output_dir)


    ##
    ## run plink r2
    ##



    # if not glob.glob(snps_dir+"/chr*"):
    #     logger.error("SNPs to run pairwise LD were not found!\nLooked in {}".format(snps_dir))
    #     sys.exit()

    keep_warnings = []

    #TO DO
    # multithread this

    # loop through each chromosome file that has one row per lead SNP to calculate r2 for
    for this_chr_num, lead_snps_file in snp_list_by_chr_files_dict.items():
    # for this_chr_num in np.arange(1,23):
    # for snps_in_this_chr in glob.glob(snps_dir+"/chr*"):

        # this_chr = os.path.split(snps_in_this_chr)[1].split("_")[0]  # e.g. chr1
        # this_chr_num = this_chr.split('chr')[1]
        #

        # set up plink
        outputfile = os.path.join(plink_output_dir, "chr{}".format(this_chr_num))
        this_chr_gwas_snp = gwas_snps_by_chr_files_dict[this_chr_num]

        # here we extract all gwas snps in this chromosome and ld expand the lead snps
        plinkcmd = f'plink --bfile {thous_gen_file.format(this_chr_num)} --extract {this_chr_gwas_snp} --ld-snp-list {lead_snps_file} --r2 inter-chr --ld-window-r2 {min_r2_threshold} --out {outputfile}'
        logger.debug(plinkcmd)

        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()


        _ = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        keep_warnings.append(any_warnings)

    if np.any(keep_warnings):
        print("There were warnings while running plink --r2")





    # parse ld file & write summary
    if not glob.glob(plink_output_dir+"/chr*.ld"):
        logger.error("Pairwise LD files were not found!\n Looked in {}".format(plink_output_dir))

    store_ld_df = pd.DataFrame()
    for thischr in glob.glob(plink_output_dir+"/chr*.ld"):

        clump_df = pd.read_csv(thischr, sep="\s+")
        store_ld_df = store_ld_df.append(clump_df.loc[:, ['SNP_A', 'SNP_B', 'R2']])


    store_ld_df['snpA_B'] = store_ld_df['SNP_A'] + "_" + store_ld_df['SNP_B']
    store_ld_df['snpB_A'] = store_ld_df['SNP_B'] + "_" + store_ld_df['SNP_A']
    store_ld_df.to_csv(ld_summary_file, index=False, header=True, sep="\t")

    logger.debug(f"Pairwise LD for all input snps saved to: {ld_summary_file}")
    logger.info("Done calculing pariwise LD. It took {0:.2f} minutes.".format((time.time() - tstart)/60))


    return store_ld_df



def calc_r2_among_snps_in_snp_list(snp_list_by_chr_files_dict, output_dir, thous_gen_file, min_r2_threshold=0.7):
    """"Calculates pariwise LD (r2) using plink --r2 inter-chr. See below for parameter details
        INPUTS:
            *
            * output_dir: base directory where all resutls will be store; if it doesn't exist, it will create.
                > subdirectories will not be overwritten, will exit script if it already exists.
            * thous_gen_file: path to Thousand Genomes PLINK file; no further filtering is done
        OUTPUTS:
            * plink output is written to: plink_ld_output/.
            * pairwise ld is combined across chromosomes and written to: all_pairwise_ld.tsv
            * plink screen output is written to: plink_ld_stdout.txt"""

    ##
    ## setup
    ##

    tstart = time.time()

    # set up output folders and files
    ld_summary_file = os.path.join(output_dir, "all_pairwise_ld.tsv")
    plink_output_dir = os.path.join(output_dir, 'plink_ld_output')


    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    safe_mkdir(plink_output_dir)


    ##
    ## run plink r2
    ##



    if len(snp_list_by_chr_files_dict) == 0:
        logger.error("SNPs to run pairwise LD were not found!")
        sys.exit()

    keep_warnings = []

    #TO DO
    # multithread this

    for this_chr_num, snps_in_this_chrom_file_path in snp_list_by_chr_files_dict.items():

        # set up plink
        outputfile = os.path.join(plink_output_dir, "chr{}".format(this_chr_num))
        plinkcmd = f'plink --bfile {thous_gen_file.format(this_chr_num)} --extract {snps_in_this_chrom_file_path} --r2 inter-chr --ld-window-r2 {min_r2_threshold} --out {outputfile}'
        logger.debug(plinkcmd)

        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()


        _ = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        keep_warnings.append(any_warnings)

    if np.any(keep_warnings):
        print("There were warnings while running plink --r2, check logs.")





    # parse ld file & write summary
    if not glob.glob(plink_output_dir+"/chr*.ld"):
        logger.error("Pairwise LD files were not found!\n Looked in {}".format(plink_output_dir))

    store_ld_df = pd.DataFrame()
    for thischr in glob.glob(plink_output_dir+"/chr*.ld"):

        clump_df = pd.read_csv(thischr, sep="\s+")
        store_ld_df = store_ld_df.append(clump_df.loc[:, ['SNP_A', 'SNP_B', 'R2']])


    if store_ld_df.shape[0] == 0:

        logger.debug(f"No LD snps were found after running plink.")
        logger.debug("Done calculing pariwise LD. It took {0:.2f} minutes.".format((time.time() - tstart)/60))

        return store_ld_df

    else:
        store_ld_df['snpA_B'] = store_ld_df['SNP_A'] + "_" + store_ld_df['SNP_B']
        store_ld_df['snpB_A'] = store_ld_df['SNP_B'] + "_" + store_ld_df['SNP_A']
        store_ld_df.to_csv(ld_summary_file, index=False, header=True, sep="\t")

        logger.debug(f"Pairwise LD for all input snps saved to: {ld_summary_file}")
        logger.debug("Done calculing pariwise LD. It took {0:.2f} minutes.".format((time.time() - tstart)/60))


    return store_ld_df

def ld_expand_snp_list(snp_list_by_chr_files_written, plink_r2_output_dir, thous_gen_file, r2_ldexp_threshold = 0.7):

    keep_warnings= []
    tstart = time.time()
    for this_chr_num, snps_in_this_chr_file in snp_list_by_chr_files_written.items():


        # set up plink to ld expand out to r2_ldexp_threshold

        plink_outputfile = os.path.join(plink_r2_output_dir, "chr{}".format(this_chr_num))
        # --ld-window-r2 is < specified threshold r2 will not be reported
        plinkcmd = f'plink --bfile {thous_gen_file.format(this_chr_num)} --ld-snp-list {snps_in_this_chr_file} --r2 inter-chr --ld-window-r2 {r2_ldexp_threshold} --out {plink_outputfile}'
        logger.debug(plinkcmd)


        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()

        any_errors = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        keep_warnings.append(any_warnings)


    # error check the plink run
    if any_errors:
        sys.exit("There were errors while running plink --r2 to ld expand input snps.")


    if np.any(keep_warnings):
        print("There were warnings while running plink --r2")



    # combine together into one file
    if not glob.glob(plink_r2_output_dir+"/chr*.ld"):
        logger.error("Pairwise LD files were not found!\n Looked in {}".format(plink_r2_output_dir))

    store_ld_df = pd.DataFrame()
    for thischr in glob.glob(plink_r2_output_dir+"/chr*.ld"):

        clump_df = pd.read_csv(thischr, sep="\s+")
        store_ld_df = store_ld_df.append(clump_df.loc[:, ['SNP_A', 'SNP_B', 'R2']])


    store_ld_df['snpA_B'] = store_ld_df['SNP_A'] + "_" + store_ld_df['SNP_B']
    store_ld_df['snpB_A'] = store_ld_df['SNP_B'] + "_" + store_ld_df['SNP_A']

    store_ld_df.reset_index(inplace=True, drop=True)
    logger.debug("Done calculing pariwise LD. It took {0:.2f} minutes.".format((time.time() - tstart)/60))

    return store_ld_df