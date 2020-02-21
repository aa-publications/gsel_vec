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


# def calc_r2_for_input_snps(snps_dir, output_dir, thous_gen_file, logger ):
def calc_r2_for_input_snps(snps_dir, output_dir, thous_gen_file ):
    """"Calculates pariwise LD (r2) using plink --r2 inter-chr. See below for parameter details

        INPUTS:
            * snps_dir: directory where there is one file per chromosome with one variant IDs per line
                > variantID should be in CHR:BP format
                > these are the input snps for which pairwise LD will be calcualted

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
    plink_screenoutput_file = os.path.join(output_dir, "plink_ld_stdout.txt")

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    safe_mkdir(plink_output_dir)


    ##
    ## run plink r2
    ##

    fh_plinkoutput = open(plink_screenoutput_file, 'w')

    if not glob.glob(snps_dir+"/chr*"):
        logger.error("SNPs to run pairwise LD were not found!\nLooked in {}".format(snps_dir))
        sys.exit()

    keep_warnings = []

    #TO DO
    # multithread this


    for snps_in_this_chr in glob.glob(snps_dir+"/chr*"):

        this_chr = os.path.split(snps_in_this_chr)[1].split("_")[0]  # e.g. chr1
        this_chr_num = this_chr.split('chr')[1]

        # set up plink
        outputfile = os.path.join(plink_output_dir, "{}".format(this_chr))
        plinkcmd = f'plink --bfile {thous_gen_file.format(this_chr_num)} --extract {snps_in_this_chr} --r2 inter-chr --ld-window-r2 0.7 --out {outputfile}'

        # TODO: delete this... if there are any spaces in file paths, the script will break...
        tempcmd = plinkcmd.replace('Google Drive', 'Google_Drive')
        plinkcmd=[x.replace('Google_Drive','Google Drive') for x in tempcmd.split()]

        logger.debug(' '.join(plinkcmd))

        # plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_process = Popen(plinkcmd, stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()

        _ = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        keep_warnings.append(any_warnings)

    if np.any(keep_warnings):
        print("There were warnings while running plink --r2")


    fh_plinkoutput.close()


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

    logger.debug(f"Pairwise LD for all GWAS input snps saved to: {ld_summary_file}")
    logger.info("Done calculing pariwise LD. It took {0:.2f} minutes.".format((time.time() - tstart)/60))


    return store_ld_df