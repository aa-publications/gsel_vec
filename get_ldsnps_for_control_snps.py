#!/bin/python
# - input is the file with control SNPs
# - runs plink to find ld snps for all control snps
# - outputs this by chromosome
#
# Abin Abraham
# created on: 'now'


import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from subprocess import Popen, PIPE

from helper_general import Outputs, error_check_plink_run, warning_check_plink_run


import logging
logger = logging.getLogger('main.{}'.format(__name__))




# PATHS, can remove...
# THOU_GENOME_DIR="/dors/capra_lab/projects/gwas_allele_age_evolution/data/1KG_snpsnap/nodups"
# THOU_GENOME_DIR = '1kg'
# thous_gen_file = os.path.join(THOU_GENOME_DIR, "EUR.chr{}.phase3.nodups")



DATE = datetime.now().strftime('%Y-%m-%d')


# -----------
# FUNCTIONS
# -----------

def set_up_outputs(OutputObj):

    # set up ouput dirs
    OutputObj.add_output('temp_ldsnps_dir', "temp_ldsnps_to_run_plink", mkdir=True, add_root=True)
    OutputObj.add_output('ld_snps_for_control_snps_dir', "ld_snps_for_control_snps", mkdir=True, add_root=True)

    # set up ouput files
    OutputObj.add_output('ldsnp_list_file', "ldsnps_chr{}.txt", custom_root=OutputObj.get('temp_ldsnps_dir'))
    OutputObj.add_output('plink_output_file', "ld_for_control_snps_chr{}", custom_root=OutputObj.get('ld_snps_for_control_snps_dir'))



    return OutputObj


def parse_input_args():
    parser = argparse.ArgumentParser(description='Given a set of control snps, this script will run plink to calculate snps in ld.')


    # TODO: delete (only for dev purposes)
    if len(sys.argv) != 5:
        csnps_file = "test_output/input_gwas_file.tsv"
        output_root= os.getcwd()
    else:

        ### REQUIRED ARGUMENTS IN ORDER
        parser.add_argument('control_file', action='store', type=str, help="output from match_snps.py")

        parser.add_argument('one_kg_file', action='store', type=str, help="full path to 1KG file with open template for chromosome")
        parser.add_argument('output_file', action='store', type=str)


        # control_file
        #  -- this is the output from match_snps.py
        # tsv file with one row per snp
        # first column is the input gwas snp
        # remaining columns are 'control sets' (one per column) having control snps

        results = parser.parse_args()

        #retrieve passed arguments
        csnps_file = results.control_file
        thous_gen_file = results.one_kg_file
        output_root = results.output_file




    return csnps_file, thous_gen_file, output_root


def get_ldsnps_for_control_snps(control_snps_file, thous_gen_file, output_root,  clump_r2_threshold=0.8):

    sstart = time.time()
    logger.info("Starting to find LD snps for control snps.")

    # set up outputs
    output_dir = os.path.join(output_root, 'ld_snps_for_control_snps')
    OutObj = Outputs(output_dir, overwrite=True)
    OutputObj = set_up_outputs(OutObj)



    # load all ld control snps
    snps_df = pd.read_csv(control_snps_file, sep="\t")

    # convert wide to long format
    wide_snp_df = pd.wide_to_long(snps_df, stubnames='Set_', i="lead_snp", j="Set")
    wide_snp_df.reset_index(inplace=True)
    wide_snp_df['chromosome']  = wide_snp_df.Set_.apply(lambda x: int(x.split(":")[0]))
    wide_snp_df.sort_values('Set_', inplace=True)

    # 3) write unique ld control snps to file
    for chrm in np.arange(1,23,1):

        snps_to_write = wide_snp_df.loc[wide_snp_df['chromosome'] ==chrm, 'Set_' ].unique().tolist()

        if len(snps_to_write) > 0:
            with open(OutputObj.get('ldsnp_list_file').format(chrm), 'w') as fw:
                for ldsnp in snps_to_write:
                    fw.write(ldsnp + "\n")


    # 4) run plink r2
    for chrm in np.arange(1,23,1):

        # skip chromosome if the file does not exists...
        if not os.path.isfile(OutputObj.get('ldsnp_list_file').format(chrm)):
            continue

        # ld-window-r2  --> minimum r2 requrieed to report in output
        plinkcmd = f"plink --bfile {thous_gen_file.format(chrm)} --r2 dprime --ld-snp-list {OutputObj.get('ldsnp_list_file').format(chrm)} --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 {clump_r2_threshold} --out {OutputObj.get('plink_output_file').format(chrm)}"


        # 5) start plink cmd
        pstart = time.time()
        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()

        error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)

        pend = time.time()
        logger.debug(">>>Done with chr{}. Took {:.2f} mins.".format( chrm, (pend-pstart)/60))


    logger.info("Done getting LD snps for control snps. It took {:.2f} minutes.".format( (time.time() - sstart)/60))

    return OutputObj



# -----------
# MAIN
# -----------

if __name__ == "__main__":

    control_snps_file, thous_gen_file,  output_root = parse_input_args()

    get_ldsnps_for_control_snps(control_snps_file, thous_gen_file, output_root)

