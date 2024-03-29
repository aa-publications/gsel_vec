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


from .helper_calc_r2 import (
    calc_r2_for_input_snps,
    calc_r2_among_snps_in_snp_list,
    ld_expand_snp_list,
)  # pylint:E0401
from .helper_general import Outputs, error_check_plink_run, warning_check_plink_run
from .helper_clump_snps import (
    write_gwas_sumstats_by_chr,
    get_list_of_ld_snps,
    get_r2_for_lead_ld_snps,
    bin_ldsnp_per_leadsnp,
    write_snp_list_by_chr,
)

tstart = time.time()

from datetime import datetime

DATE = datetime.now().strftime("%Y-%m-%d")


import logging

logger = logging.getLogger("main.{}".format(__name__))

# %%


# -----------
# DEPENDENCIES
# -----------


# -----------
# FUNCION
# -----------


def parse_log_for_err(plink_clump_output_dir, missing_snps_file, plink_log_errors_file):

    store_missing_variants = list()
    store_errors = list()

    for thislog in glob.glob(plink_clump_output_dir + "/chr*.log"):
        with open(thislog, "r") as fh_log:

            for line in fh_log:
                if line.startswith("Warning:") and line.endswith("top variant.\n"):
                    missing_variant = line.split()[1][1:-1]
                    store_missing_variants.append(missing_variant)

                if line.startswith("Error:"):
                    store_errors.append("{}: {}".format(thislog, line))

    if store_missing_variants:
        with open(missing_snps_file, "w") as fout:
            fout.write("\n".join(store_missing_variants))

    if store_errors:
        with open(plink_log_errors_file, "w") as ferr:
            ferr.write("\n".join(store_errors))

    return store_errors, store_missing_variants


def concat_plink_clump_output(plink_clump_output_dir):

    store_clump_df = pd.DataFrame()
    for thisclump in glob.glob(plink_clump_output_dir + "/chr*.clumped"):

        clump_df = pd.read_csv(thisclump, sep="\s+")
        store_clump_df = store_clump_df.append(clump_df)

    if store_clump_df.shape[0] == 0:
        raise RuntimeError(f"Plink clumped outputs file are empty...")

    # format and write concatenated clump file
    store_clump_df["CHR_BP"] = (
        store_clump_df["CHR"].map(str) + ":" + store_clump_df["BP"].map(str)
    )

    return store_clump_df


def parse_input_args():

    # TODO: delete (only for dev purposes)
    if len(sys.argv) != 5:
        print("running dev args...delete later...")
        gwas_summary_file = "test_data/input_gwas_file.tsv"
        output_root = os.getcwd()
    else:

        parser = argparse.ArgumentParser(
            description="Will LD clump from GWAS summary stats."
        )

        parser.add_argument(
            "-g",
            "--gwas_summary_stat_file",
            dest="gwas_file",
            action="store",
            type=str,
            help="path to gwas summary file",
        )

        parser.add_argument(
            "-k",
            "--1kg_file",
            dest="one_kg_file",
            action="store",
            type=str,
            help="path to gwas summary file",
        )

        parser.add_argument(
            "-o",
            "--output",
            dest="output_root",
            action="store",
            type=str,
            help="path to dir where output will be stored",
        )

        argparser = parser.parse_args()
        gwas_summary_file = argparser.gwas_file
        thous_gen_file = argparser.one_kg_file
        output_root = argparser.output_root

    return gwas_summary_file, thous_gen_file, output_root


def set_up_outputs(OutputObj):

    # set up ouput dirs
    OutputObj.add_output("by_chr_dir", "gwas_stats_by_chr", mkdir=True, add_root=True)
    OutputObj.add_output(
        "gstats_pos_by_chr_dir", "gwas_stats_pos_only_by_chr", mkdir=True, add_root=True
    )
    OutputObj.add_output(
        "plink_clump_output_dir", "plink_clump_output", mkdir=True, add_root=True
    )
    OutputObj.add_output(
        "gstats_r2_dir", "pairwise_r2_between_gwas_snps", mkdir=True, add_root=True
    )
    OutputObj.add_output(
        "lead_snps_pos_by_chr_dir",
        "lead_snps_pos_by_chr_dir",
        mkdir=True,
        add_root=True,
    )
    # set up ouput files

    OutputObj.add_output("reprint_gwas_file", "input_gwas_file.tsv", add_root=True)
    OutputObj.add_output(
        "missing_snps_file", "input_gwas_snps_not_found.txt", add_root=True
    )
    OutputObj.add_output("lead_snps_file", "lead_gwas_snps.tsv", add_root=True)
    OutputObj.add_output(
        "plink_log_errors_file", "plink_clump_errors.txt", add_root=True
    )
    OutputObj.add_output("clumped_snps_file", "gwas_snps_clumped.tsv", add_root=True)
    OutputObj.add_output(
        "bin_by_ld_file", "lead_gwas_snps_with_ldsnps_counts.txt", add_root=True
    )
    OutputObj.add_output(
        "lead_snp_ld_pairs_r2", "lead_gwas_snp_and_ld_snp_r2.tsv", add_root=True
    )
    OutputObj.add_output(
        "log_file", "log_clump_snps_{}.txt".format(DATE), add_root=True
    )

    return OutputObj


def set_up_outputs_for_snp_list(OutObj):
    OutObj.add_output("by_chr_dir", "snps_by_chr", mkdir=True, add_root=True)

    # OutObj.add_output('ld_expand_r2_dir', "ld_expand_snp_list", mkdir=True, add_root=True)
    OutObj.add_output(
        "plink_r2_output", "ld_expand_plink_output", mkdir=True, add_root=True
    )
    OutObj.add_output(
        "all_input_snps_ld_expanded", "all_input_snps_ld_expanded.tsv", add_root=True
    )

    OutObj.add_output(
        "pairwise_r2_of_input_snps",
        "pairwise_r2_of_input_snps",
        add_root=True,
        mkdir=True,
    )
    OutObj.add_output(
        "bin_by_ld_file", "input_snps_with_binned_ldsnps_counts.tsv", add_root=True
    )
    OutObj.add_output("lead_snp_ld_pairs_r2", "input_and_ld_snps_r2.tsv", add_root=True)
    OutObj.add_output("non_indep_snps_file", "non_indep_input_snps.tsv", add_root=True)

    # set up ouput files
    OutObj.add_output(
        "missing_snps_file", "input_snps_not_in_1kg_excluded.txt", add_root=True
    )
    OutObj.add_output(
        "plink_log_errors_file", "plink_ld_expand_errors.txt", add_root=True
    )

    return OutObj


def force_input_snp_in_first_col(keep_autosomal_snps, og_store_ld_df):
    """ ensure that the input snp is in SNP_A column; if there are no LD snps, then there is a row for that SNP with 'NONE' as placeholders """

    store_ld_df = og_store_ld_df.copy()

    # check that all input snps are in SNP_A column
    remaining_autosomal_snps_in_SNP_A = set(keep_autosomal_snps).difference(
        set(store_ld_df.SNP_A.unique())
    )
    input_autosomal_snps_in_SNP_B = (
        set(keep_autosomal_snps)
        .difference(set(store_ld_df.SNP_B.unique()))
        .difference(remaining_autosomal_snps_in_SNP_A)
    )  # remove any snps found in SNP_A column already

    # if  input snps that are not in SNP_A
    if len(remaining_autosomal_snps_in_SNP_A) != 0:

        # check if any missing input snps from SNP_A are in SNP_B column
        remaining_input_snps_in_snp_b = set(
            remaining_autosomal_snps_in_SNP_A
        ).difference(set(store_ld_df.SNP_B.unique()))

        if len(remaining_input_snps_in_snp_b) > 0:

            # SWAP THE SNP_A and SNP_B labels
            temp_snpb = store_ld_df.loc[
                store_ld_df["SNP_B"].isin(remaining_input_snps_in_snp_b)
            ].copy()
            assert temp_snpb.columns.values.tolist() == [
                "SNP_A",
                "SNP_B",
                "R2",
                "snpA_B",
                "snpB_A",
            ], "columns are not what is expected"
            temp_snpb.columns = ["SNP_B", "SNP_A", "R2", "snpB_A", "snpA_B"]

            # drop these snps from the initial df and add the swapped columsn back in
            store_ld_df.drop(temp_snpb.index, axis=0, inplace=True)

            # add the swapped snps back in
            store_ld_df = store_ld_df.append(temp_snpb, sort=True)

        # check again if there are missing snps from SNP_A column
        remaining_autosomal_snps_in_SNP_A_after_swap = set(
            keep_autosomal_snps
        ).difference(set(store_ld_df.SNP_A.unique()))

        if len(remaining_autosomal_snps_in_SNP_A_after_swap) > 0:

            add_miss_snp_df = pd.DataFrame(
                {
                    "SNP_A": list(remaining_autosomal_snps_in_SNP_A_after_swap),
                    "SNP_B": "NONE",
                    "R2": "NONE",
                    "snpA_B": "NONE",
                    "snpB_A": "NONE",
                }
            )
            store_ld_df = store_ld_df.append(add_miss_snp_df, sort=True)

    # reorder columns
    store_ld_df = store_ld_df.loc[:, ["SNP_A", "SNP_B", "snpA_B", "snpB_A", "R2"]]
    return store_ld_df


def clump_snps(
    gwas_summary_file,
    output_root,
    thous_gen_file,
    lead_snp_min_gwas_pvalue=0.00000005,
    ld_snps_min_gwas_pvalue=0.00000005,
    min_r2_to_clump=0.9,
    min_kb_from_index_to_clump=250,
    ld_expand_lead_snp_min_r2=0.7,
):

    output_dir = os.path.join(output_root, "ld_clump_inputs")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)

    # TODO: figure out how to deal with input data...
    # excpected header names for gwas summary statistics file
    # ALSO: we expect the chromosome column to only have the chromsome number e.g. 2, 3, 4, not CHR2
    gsum_header = {"rsID": "snp", "chr": "chr", "basepair": "pos", "pvalue": "p"}

    ##
    ## parse gwas input snps
    ##

    raw_gwas_df = pd.read_csv(
        gwas_summary_file,
        sep="\t",
        dtype={"snp": str, "chr": np.int, "pos": np.int, "p": np.float64},
    )

    raw_gwas_df.to_csv(
        OutObj.get("reprint_gwas_file"), sep="\t", index=False
    )  # write a copy of the input data
    logger.info(f"[status] Loaded GWAS summary stats with {raw_gwas_df.shape[0]} rows.")

    # write input snps by chromosome
    logger.debug(f"Splitting GWAS summary stats by chromosome.")
    gwasfilename = os.path.split(gwas_summary_file)[1]

    ###
    ###    error check gwas summary stats
    ###

    # check for duplicate positions based on chr and pos
    dup_rows = raw_gwas_df.duplicated(subset=["chr", "pos"], keep="first")
    gwas_df = raw_gwas_df.loc[~dup_rows].copy()
    logger.info(
        f"{np.sum(dup_rows)} duplciated variants (based on chr and pos) were removed"
    ) if np.any(dup_rows) else None

    # remove rows whose chromosome is not 1 through 22
    gwas_df = gwas_df.loc[gwas_df["chr"].isin(np.arange(1, 23)), :].copy()

    ##
    ## run plink clump by chromosome
    ##
    (
        gwas_stats_by_chr_files_dict,
        gwas_snps_by_chr_files_dict,
    ) = write_gwas_sumstats_by_chr(
        gwas_df,
        OutObj.get("by_chr_dir"),
        OutObj.get("gstats_pos_by_chr_dir"),
        gsum_header["chr"],
        gsum_header["basepair"],
        gwasfilename,
        gsum_header["pvalue"],
    )

    logger.debug(f"Running plink clump on gwas variants by chromosome.")
    keep_warnings = []
    # for thisfile in glob.glob(OutObj.get('by_chr_dir')+"/chr*"):
    for chr_num, thisfile in gwas_stats_by_chr_files_dict.items():

        # inputs and outputs for plink call
        # this_chr = os.path.split(thisfile)[1].split("_")[0]  # e.g. "chr3"
        # chr_num = this_chr.split('chr')[1]
        outputfile = os.path.join(
            OutObj.get("plink_clump_output_dir"), "chr{}".format(chr_num)
        )

        plinkcmd = f"plink --bfile {thous_gen_file.format(chr_num)} --clump {thisfile} --clump-kb {min_kb_from_index_to_clump} --clump-r2 {min_r2_to_clump} --clump-p1 {lead_snp_min_gwas_pvalue} --clump-p2 {ld_snps_min_gwas_pvalue} --out {outputfile}"
        logger.debug(plinkcmd)

        # run plink cmd
        plink_process = Popen(plinkcmd.split(), stdout=PIPE, stderr=PIPE)
        plink_stdout, plink_stderr = plink_process.communicate()

        _ = error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=logger)
        any_warnings = warning_check_plink_run(
            plink_stdout, plink_stderr, plinkcmd, logger=logger
        )
        keep_warnings.append(any_warnings)


    #TODO: add a check to see if """ Warning: No significant --clump results.  Skipping"""

    if np.any(keep_warnings):
        logger.debug("Plink threw warnings...check the log file.")

    if len(glob.glob(OutObj.get("plink_clump_output_dir") + "/*.clumped")) == 0:
        raise RuntimeError(
            f'No plink clumped results in:\n{OutObj.get("plink_clump_output_dir")}'
        )
    else:
        _, store_missing_variants = parse_log_for_err(
            OutObj.get("plink_clump_output_dir"),
            OutObj.get("missing_snps_file"),
            OutObj.get("plink_log_errors_file"),
        )

    # write all plink clump output by chromosome to one file
    store_clump_df = concat_plink_clump_output(OutObj.get("plink_clump_output_dir"))
    store_clump_df.rename(columns={"SNP": "lead_snp"}, inplace=True)
    store_clump_df.to_csv(
        OutObj.get("clumped_snps_file"), index=False, header=True, sep="\t"
    )


    logger.debug(
        "[status] Plink clump found {} lead snps.".format(store_clump_df.shape[0])
    )

    ##
    ## calc r2 between gwas lead and LD SNPs
    ##

    logger.debug(
        "Calculating pairwise r2 for lead and LD gwas variants based on clumping parameters."
    )

    # load all lead_snp and ld_snp pairs
    lead_ld_df = pd.DataFrame()
    lead_snps_per_chr_files_dict = dict()
    for chrom_num in np.arange(1, 23):

        clumped_file = os.path.join(
            OutObj.get("plink_clump_output_dir"), "chr{}.clumped".format(chrom_num)
        )

        if os.path.isfile(clumped_file):

            # get all LD snps with lead SNPs (the ld snps are determined by clumping parameters specficied earlier..)
            # lead snps w/o LD snps will be included
            temp_df = pd.DataFrame(
                get_list_of_ld_snps(clumped_file), columns=["lead_snp", "ld_snp"]
            )
            lead_ld_df = lead_ld_df.append(temp_df)

            # write lead snps per chromsome
            pos_out_file = os.path.join(
                OutObj.get("lead_snps_pos_by_chr_dir"),
                os.path.basename(clumped_file).split(".clumped")[0] + ".lead_snps",
            )
            lead_snps_per_chr_files_dict[chrom_num] = pos_out_file

            lead_snps_to_write = temp_df.lead_snp.unique()
            with open(pos_out_file, "w") as fw:
                for lead_snp_to_write in lead_snps_to_write:
                    fw.writelines(lead_snp_to_write + "\n")

    # check if there are no LD snps ...
    if lead_ld_df.shape[0] == 0:
        logger.debug(
            "[status] No lead or ld snps found after clumping."
        )
        sys.exit("No lead or ld snps found after clumping.")

    # for lead snps, LD expand them
    # note: if lead SNP doesn't have any LD snps at specified r2, then it will not be reported in plink's outputs ( take care of this downstream in the pipeline)
    # note: min_r2_threshold is the r2 to which the lead snp should be expanded to (this just has tobe a lower or equal r2 threshold than what will be used to ld expand the control snps)
    r2_df = calc_r2_for_input_snps(
        lead_snps_per_chr_files_dict,
        gwas_snps_by_chr_files_dict,
        OutObj.get("gstats_r2_dir"),
        thous_gen_file,
        min_r2_threshold=ld_expand_lead_snp_min_r2,
    )

    ###
    ###    bin ld snps by LD deciles
    ###

    # get r2 for all lead and ld snp pairs
    #   note: lead snps w/o any LD snps in r2_df will also have a row with NONE as a placeholder

    store_lead_snp_ld_df = get_r2_for_lead_ld_snps(lead_ld_df, r2_df)
    store_ld_bins_df = bin_ldsnp_per_leadsnp(store_lead_snp_ld_df)
    store_ld_bins_df.reset_index(drop=False, inplace=True)
    store_ld_bins_df.rename(columns={"index": "lead_snp"}, inplace=True)

    ##
    ## write outputs
    ##

    # lead snps only file
    # TODO: delete the label to this file from the OUTOBJ object
    # lead_ld_df[~lead_ld_df.duplicated('lead_snp', keep='first')].to_csv(OutObj.get('lead_snps_file'), sep="\t", index=False, header=False, columns=['lead_snp'])
    # logger.info(f"Wrote lead snps only: {OutObj.get('lead_snps_file')}")

    # ld_bins
    store_ld_bins_df.to_csv(OutObj.get("bin_by_ld_file"), sep="\t", index=False)
    logger.debug(
        f"Wrote table of number of LD snps (binned): {OutObj.get('bin_by_ld_file')}"
    )

    # lead snp, LD snp, r2 table
    store_lead_snp_ld_df.to_csv(
        OutObj.get("lead_snp_ld_pairs_r2"), sep="\t", index=False
    )
    logger.debug(
        f"Wrote table of lead and ld snp pairs with r2: {OutObj.get('lead_snp_ld_pairs_r2')}"
    )

    logger.info(
        f"* While running plink clump, {len(store_missing_variants)} input SNPs not found in 1KG written."
    )
    logger.debug(
        f"[status] Done clumping GWAS summary stats. Found {store_ld_bins_df.lead_snp.nunique() } lead snps and took {(time.time()-tstart)/60:.2f} minutes."
    )

    return OutObj


def clump_snp_list(
    snps_list_file,
    output_root,
    thous_gen_file,
    min_r2_to_ld_exp_lead_snp=0.9,
    min_r2_for_input_snp_indep=0.9,
):

    ostart = time.time()
    # set up outputs
    output_dir = os.path.join(output_root, "ld_clump_inputs")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs_for_snp_list(OutObj)

    ###
    ###    load list of snp and check input format
    ###

    # load input snps
    with open(snps_list_file, "r") as fo:
        snp_file_contents = fo.read().splitlines()

    logger.info(f"Loaded {len(snp_file_contents):,} snps for analysis.")

    # ensure input snps follow the following format chr:pos
    assert np.all(
        [
            True if re.match(r"[1-9][0-9]*:[0-9]+", str(snp)) else False
            for snp in snp_file_contents
        ]
    ), "check input snps are in formatted correctly (2:22222, chr:pos, without any prefix for chromosome)"

    # keep only autosomal snps
    keep_autosomal_snps = [
        snp for snp in snp_file_contents if (int(snp.split(":")[0]) < 23)
    ]
    non_auto_snps_removed = set(snp_file_contents).difference(set(keep_autosomal_snps))

    if len(non_auto_snps_removed) > 0:
        logger.info(
            f"Removed {len(non_auto_snps_removed)} snps for not being autosomal snps."
        )

    if len(keep_autosomal_snps) == 0:
        sys.exit(
            "Error, script terminated because no autosomal snps left in the input snp list."
        )

    ###
    ###    ld expand each input snp
    ###

    # split by chromosome
    snp_list_by_chr_files_dict = write_snp_list_by_chr(
        keep_autosomal_snps, OutObj.get("by_chr_dir")
    )

    # ld expand
    plink_r2_output_dir = OutObj.get("plink_r2_output")
    all_snps_ld_expanded_file = OutObj.get("all_input_snps_ld_expanded")
    store_ld_df = ld_expand_snp_list(
        snp_list_by_chr_files_dict,
        plink_r2_output_dir,
        thous_gen_file,
        r2_ldexp_threshold=min_r2_to_ld_exp_lead_snp,
    )  # r2_ldexp_threshold is < threshold will not be reported

    _, store_miss_variants = parse_log_for_err(
        plink_r2_output_dir,
        OutObj.get("missing_snps_file"),
        OutObj.get("plink_log_errors_file"),
    )

    ###
    ###    bin snps by ld
    ###

    # remove snps that were not found in 1kg database (so that these SNPs don't get filled in as having no SNP in LD with it downstream)
    keep_autosomal_snps = set(keep_autosomal_snps).difference(set(store_miss_variants))

    if len(keep_autosomal_snps) == 0:
        sys.exit(
            "Error, script terminated because no autosomal snps left in the input snp list."
        )
    else:
        logger.info(
            f"{len(keep_autosomal_snps):,} input snps remains after intersecting with thousand Genomes PhaseIII."
        )

    # remove rows where SNP_A == SNP_B and R2 ==1
    store_ld_df.drop(
        store_ld_df[
            (store_ld_df["SNP_A"] == store_ld_df["SNP_B"]) & (store_ld_df["R2"] == 1)
        ].index,
        axis=0,
        inplace=True,
    )

    # force input snp in SNP_A and bin, makes sure that each snp in keep_autosomal_snps have at least one row
    for_ld_binning_df = force_input_snp_in_first_col(keep_autosomal_snps, store_ld_df)
    for_ld_binning_df.rename(columns={"SNP_A": "lead_snp"}, inplace=True)
    store_ld_bins_df = bin_ldsnp_per_leadsnp(for_ld_binning_df)
    store_ld_bins_df.reset_index(drop=False, inplace=True)
    store_ld_bins_df.rename(columns={"index": "lead_snp"}, inplace=True)

    ###
    ###    calculate ld between input snps and check for ld independence
    ###
    r2start = time.time()
    pairwise_r2_df = calc_r2_among_snps_in_snp_list(
        snp_list_by_chr_files_dict,
        OutObj.get("pairwise_r2_of_input_snps"),
        thous_gen_file,
        min_r2_threshold=min_r2_for_input_snp_indep,
    )
    logger.info(
        f"Done calculating r2 between input snps. Took {(time.time()-r2start)/60:.2f} minutes."
    )

    # write a list of input snps are within the given r2 threshold
    if pairwise_r2_df.shape[0] > 0:
        forced_pairwise_r2_df = force_input_snp_in_first_col(
            keep_autosomal_snps, pairwise_r2_df
        )
        non_indep_input_snps_df = pairwise_r2_df.loc[
            pairwise_r2_df["R2"] > min_r2_for_input_snp_indep
        ].copy()
        write_non_indep_df = non_indep_input_snps_df.loc[
            non_indep_input_snps_df["SNP_A"].isin(keep_autosomal_snps)
        ].copy()
        write_non_indep_df.to_csv(
            OutObj.get("non_indep_snps_file"), sep="\t", index=False, header=True
        )
        logger.debug(
            f"Wrote input snps that are non-indepdent at r2 > {min_r2_for_input_snp_indep} to: {OutObj.get('non_indep_snps_file')}"
        )

    ###
    ###    write
    ###

    # write input snps and ld snps
    #   - seperate lead snps with and without ≥1 ld snp
    lead_snp_w_no_ld_snps_df = for_ld_binning_df.loc[
        for_ld_binning_df["SNP_B"] == "NONE"
    ].copy()
    lead_snp_w_ld_snps_df = for_ld_binning_df.loc[
        for_ld_binning_df["SNP_B"] != "NONE"
    ].copy()

    # filter on r2 adn then concat back lead snps w/o LD snps
    lead_and_ld_snp_df = lead_snp_w_ld_snps_df.loc[
        lead_snp_w_ld_snps_df["R2"] > min_r2_to_ld_exp_lead_snp,
        ["lead_snp", "SNP_B", "R2"],
    ].copy()
    lead_and_no_ld_snp_df = lead_snp_w_no_ld_snps_df.loc[
        :, ["lead_snp", "SNP_B", "R2"]
    ].copy()
    all_lead_ld_snp_df = pd.concat([lead_and_ld_snp_df, lead_and_no_ld_snp_df], axis=0)

    # write
    all_lead_ld_snp_df.rename(columns={"SNP_B": "ld_snp"}, inplace=True)
    all_lead_ld_snp_df.to_csv(
        OutObj.get("lead_snp_ld_pairs_r2"), sep="\t", index=False, header=True
    )
    logger.debug(f"Wrote input and ld snp to: {OutObj.get('lead_snp_ld_pairs_r2')}")

    # write snps binned ld
    store_ld_bins_df.to_csv(OutObj.get("bin_by_ld_file"), sep="\t", index=False)
    logger.debug(
        f"Wrote table of number of LD snps (binned): {OutObj.get('bin_by_ld_file')}"
    )

    # write input snps excluded from analysis
    if store_miss_variants:
        with open(OutObj.get("missing_snps_file"), "w") as fo:
            fo.write("\n".join(store_miss_variants))

    logger.info(
        f"Done clumping input snps. Took {(time.time()-ostart)/60:.2f} minutes."
    )

    return OutObj


# -----------
# MAIN
# -----------

if __name__ == "__main__":

    gwas_summary_file, thous_gen_file, output_dir = parse_input_args()

    _ = clump_snps(gwas_summary_file, output_dir, thous_gen_file)
