#!/bin/python
# This script will mimics SNPSNAP for matching a set of input snps.
#
#
#
# Abin Abraham
# created on: 2019-12-31 12:08:19


import os
import sys
import time
import argparse

import numpy as np
import pandas as pd
from datetime import datetime


from functools import partial
from multiprocessing import Pool, cpu_count, current_process
from .helper_general import Outputs, report_mem


import logging

logger = logging.getLogger("main.{}".format(__name__))

DATE = datetime.now().strftime("%Y-%m-%d")


# REQUIRED PATHS
# SNP_ANNO_FILE="snpsnap_database/ld0.1_collection.tab.gz"
# SNP_ANNO_FILE="snpsnap_database/ld0.1_collection.tab.pickle"


# -----------
# FUNCTION
# -----------


def parse_input_args():

    parser = argparse.ArgumentParser(description="Script mimicing SNPSNAP")

    ### REQUIRED ARGUMENTS IN ORDER
    parser.add_argument(
        "input_snps_file",
        action="store",
        type=str,
        help="full path of file with one line per input snp formatted chr:pos (e.g 10:12345)",
    )

    parser.add_argument(
        "n_matches", action="store", type=int, help="number of matching snps requested"
    )
    parser.add_argument(
        "ldbud_r2",
        action="store",
        type=str,
        help="r2 threshold for LDbuddies. must be formatted as 'friends_ld05' or 'friends_ld06' etc. ",
    )

    parser.add_argument(
        "db_file",
        action="store",
        type=str,
        help="full path to snpsnap downloaded database file",
    )

    parser.add_argument(
        "output_root", action="store", type=str, help="output directory"
    )

    # retrieve passed arguments
    results = parser.parse_args()
    input_snps_file = results.input_snps_file
    n_matches = results.n_matches
    ld_buddies_r2 = results.ldbud_r2
    db_file = results.db_file
    output_root = results.output_root

    return input_snps_file, n_matches, ld_buddies_r2, db_file, output_root


def get_snps_to_match(input_snps_file):

    # getting the list of lead GWAS snps from lead and LD SNP pairs.
    # should include a row for lead SNPs without ld snps (R2 column should be "NONE")
    df = pd.read_csv(input_snps_file, sep="\t")

    snps_to_match = df.lead_snp.unique().tolist()

    return snps_to_match


def set_up_outputs(OutputObj):

    # set up ouput files
    OutputObj.add_output("excluded_snps_file", "input_snps_excluded.txt", add_root=True)
    OutputObj.add_output("matched_snps_file", "matched_snps.tsv", add_root=True)
    OutputObj.add_output(
        "quality_score_file", "match_quality_per_snp.tsv", add_root=True
    )
    OutputObj.add_output("match_summary_file", "matching_summary.txt", add_root=True)

    return OutputObj


def write_header_to_summary_file(
    fsummary, input_snps_file, n_matches, ld_buddies_r2, output_dir
):

    fsummary.writelines("Running:  \n\t{} on {}\n".format(sys.argv[0], datetime.now()))
    fsummary.writelines("Parameters:\n")
    fsummary.writelines(f"\tinput_snps_file:  {input_snps_file}\n")
    fsummary.writelines(f"\tnum_matches_requested:  {n_matches}\n")
    fsummary.writelines(f"\tr2_for_ld_buddies:  {ld_buddies_r2}\n")
    fsummary.writelines(f"\toutput_dir:  {output_dir}\n")


def get_properties_of_index_snp(snp_to_match, anno_df):
    anno_cols = ["snp_maf", "gene_count", "dist_nearest_gene", "ld_buddies"]
    index_maf, index_gene_count, index_dist_nearest_gene, index_friends_ld = (
        anno_df.loc[anno_df["snpID"] == snp_to_match, anno_cols].values[0].tolist()
    )

    return index_maf, index_gene_count, index_dist_nearest_gene, index_friends_ld


def set_thresholds(
    max_dev_maf=0.05,
    max_dev_gene_count=0.50,
    max_dev_dist_nearest_gene=0.50,
    max_dev_friends_ld=0.50,
):

    # column for each property
    # snp_maf
    # gene_count - gene density
    # dist_nearest_gene  - distance ot nearest gene
    # friends_ld08  (default)  - column based on LD

    # create 5 intervals out of max allowed deviation
    intervals_maf = np.arange(
        max_dev_maf / 5, max_dev_maf + max_dev_maf / 5, max_dev_maf / 5
    )
    intervals_percent_gene_count = np.arange(
        max_dev_gene_count / 5,
        max_dev_gene_count + max_dev_gene_count / 5,
        max_dev_gene_count / 5,
    )
    intervals_percent_dist_nearest_gene = np.arange(
        max_dev_dist_nearest_gene / 5,
        max_dev_dist_nearest_gene + max_dev_dist_nearest_gene / 5,
        max_dev_dist_nearest_gene / 5,
    )
    intervals_percent_friends_ld = np.arange(
        max_dev_friends_ld / 5,
        max_dev_friends_ld + max_dev_friends_ld / 5,
        max_dev_friends_ld / 5,
    )

    return [
        intervals_maf,
        intervals_percent_gene_count,
        intervals_percent_dist_nearest_gene,
        intervals_percent_friends_ld,
    ]


def get_thresholds(
    index_maf,
    intervals_maf,
    index_gene_count,
    intervals_percent_gene_count,
    index_dist_nearest_gene,
    intervals_percent_dist_nearest_gene,
    index_friends_ld,
    intervals_percent_friends_ld,
    attempt,
):

    maf_thresh = (
        index_maf - intervals_maf[attempt],
        index_maf + intervals_maf[attempt],
    )

    gene_count_thresh = (
        np.floor(
            index_gene_count - index_gene_count * intervals_percent_gene_count[attempt]
        ),
        np.ceil(
            index_gene_count + index_gene_count * intervals_percent_gene_count[attempt]
        ),
    )

    dist_nearest_gene_thresh = (
        np.floor(
            index_dist_nearest_gene
            - index_dist_nearest_gene * intervals_percent_dist_nearest_gene[attempt]
        ),
        np.ceil(
            index_dist_nearest_gene
            + index_dist_nearest_gene * intervals_percent_dist_nearest_gene[attempt]
        ),
    )

    friends_ld_thresh = (
        np.floor(
            index_friends_ld - index_friends_ld * intervals_percent_friends_ld[attempt]
        ),
        np.ceil(
            index_friends_ld + index_friends_ld * intervals_percent_friends_ld[attempt]
        ),
    )

    return maf_thresh, gene_count_thresh, dist_nearest_gene_thresh, friends_ld_thresh


def get_matched_snps(snp_to_match, n_matches, anno_df, thresholds):
    start = time.time()
    index_maf, index_gene_count, index_dist_nearest_gene, index_friends_ld = get_properties_of_index_snp(
        snp_to_match, anno_df
    )

    try:

        logger.debug(f"on {snp_to_match}." + report_mem())
        for attempt in range(5):

            # get interval thresholds
            intervals_maf, intervals_percent_gene_count, intervals_percent_dist_nearest_gene, intervals_percent_friends_ld = (
                thresholds
            )

            # get thresholds to match on
            maf_thresh, gene_count_thresh, dist_nearest_gene_thresh, friends_ld_thresh = get_thresholds(
                index_maf,
                intervals_maf,
                index_gene_count,
                intervals_percent_gene_count,
                index_dist_nearest_gene,
                intervals_percent_dist_nearest_gene,
                index_friends_ld,
                intervals_percent_friends_ld,
                attempt,
            )

            # get all snps in database (excluding snp of interest matching these criteria)
            matched_df = anno_df.query(
                "snp_maf >= @maf_thresh[0] and  snp_maf <= @maf_thresh[1] and gene_count >= @gene_count_thresh[0] and gene_count <= @gene_count_thresh[1] and dist_nearest_gene >= @dist_nearest_gene_thresh[0] and dist_nearest_gene <= @dist_nearest_gene_thresh[1] and ld_buddies >= @friends_ld_thresh[0] and ld_buddies <= @friends_ld_thresh[1] and snpID != @snp_to_match"
            )

            # IF THERE ARE NO MATCHES ...
            if matched_df.shape[0] == 0:
                matched_snps = ["None"] * n_matches
                num_matched_snps = 0
                return snp_to_match, matched_snps, num_matched_snps

            # randomly sample snps meeting criteria
            if matched_df.shape[0] >= n_matches:
                matched_snps = (
                    matched_df.sample(n_matches, replace=False, random_state=12)
                    .loc[:, "snpID"]
                    .values.tolist()
                )
                num_matched_snps = matched_df.shape[0]
                break

            elif attempt == 4:
                matched_snps = (
                    matched_df.sample(n_matches, replace=True, random_state=12)
                    .loc[:, "snpID"]
                    .values.tolist()
                )
                num_matched_snps = matched_df.shape[0]

    except ValueError as e:
        print(f"Value Errors in get_matched_snps while on {snp_to_match}")

    logger.info(
        "\tmatching {}. PID {}.Took {:.2f} minutes.".format(
            snp_to_match, current_process().pid, (time.time() - start) / 60
        )
    )

    return snp_to_match, matched_snps, num_matched_snps


def write_excluded_snps(snps_to_match, keep_anno_df, excluded_snps_file):

    # check if input snps are in HLA region:
    hla_snps_to_exclude = set(
        [
            x
            for x in snps_to_match
            if (int(x.split(":")[0]) == 6)
            and (25000000 < int(x.split(":")[1]) < 35000000)
        ]
    )

    # remove hla snps from input snps to match
    hla_removed_snps_to_match = set(snps_to_match).difference(hla_snps_to_exclude)

    # include snps found in the SNPSNAP database, exclude others
    snps_included = set(
        keep_anno_df.loc[
            keep_anno_df.snpID.isin(hla_removed_snps_to_match), "snpID"
        ].values.tolist()
    )
    snps_excluded = set(snps_to_match).difference(
        snps_included
    )  # will incl hla snps that were excluded

    # write excluded snps
    with open(excluded_snps_file, "w") as f_ex:
        for snp in snps_excluded:
            f_ex.write("{}\n".format(snp))

    return snps_included, snps_excluded


def write_excluded_snps_w_insuff_matches(excluded_snps_file, snps_excluded):

    # write excluded snps
    with open(excluded_snps_file, "a") as f_ex:
        for snp in snps_excluded:
            f_ex.write("{}\n".format(snp))


def write_match_quality(
    snps_to_match_ordered,
    matched_snps,
    n_matches,
    num_matched_found,
    quality_score_file,
    fsummary,
):

    # Insufficient-matches:
    #       - % of input SNPs for which SNPsnap is not able to identify the required number of matched SNPs
    # Median Match Size:
    #       - % number of matching snps found for the subset of input SNPs w/o required number of matched SNPS

    match_qual_df = pd.DataFrame(
        {
            "lead_snp": snps_to_match_ordered,
            "num_uniq_matched": num_matched_found,
            "num_match_requested": [n_matches] * len(num_matched_found),
        }
    )
    match_qual_df["insuff_matches"] = match_qual_df.num_uniq_matched < n_matches

    # mark as True if there is only < 90% of the requested number of matched snps
    # note: this is unique matched snps...
    n_match_thresh = n_matches * 0.90
    match_qual_df["n_matches_below_0.9"] = (
        match_qual_df.num_uniq_matched < n_match_thresh
    )

    match_qual_df.to_csv(quality_score_file, sep="\t", index=False)
    print("Wrote match quality per snp to: {}".format(quality_score_file))

    # calc insuff and median match proportion
    insuff_matches_percent = (
        match_qual_df.insuff_matches.sum() / match_qual_df.shape[0]
    ) * 100
    if insuff_matches_percent == 0:
        median_match_size_percent = 0
    else:
        median_match_size_percent = (
            match_qual_df.loc[match_qual_df["insuff_matches"]].num_uniq_matched.median()
            / n_matches
        ) * 100

    fsummary.writelines(
        ">>> Insufficient Match Percentage: {}\n".format(insuff_matches_percent)
    )
    fsummary.writelines(
        ">>> Median Match Size Percentage: {}\n".format(median_match_size_percent)
    )

    return match_qual_df


def match_snps(input_snps_file, n_matches, ld_buddies_r2, db_file, output_root):
    start = time.time()
    logger.info("Matching snps...")
    output_dir = os.path.join(output_root, "matched_snps_for_input_snps")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)

    # log script details
    fsummary = open(OutObj.get("match_summary_file"), "w")
    write_header_to_summary_file(
        fsummary, input_snps_file, n_matches, ld_buddies_r2, OutObj.root_dir
    )

    ###
    ###   load annotation data & format
    ###

    # load downloaded snpsnap database
    # anno_df = pd.read_pickle(db_file)
    anno_df = pd.read_csv(
        db_file,
        sep="\t",
        usecols=["snpID", "snp_maf", "gene_count", "dist_nearest_gene", ld_buddies_r2],
    )

    anno_df["ld_buddies"] = anno_df[ld_buddies_r2]  # user select r2 threshold
    smaller_anno_df = anno_df.loc[
        :, ["snpID", "snp_maf", "gene_count", "dist_nearest_gene", "ld_buddies"]
    ].copy()
    del anno_df

    # reduce memory footprint and ensure varialbe types
    smaller_anno_df.snp_maf = pd.to_numeric(
        smaller_anno_df.snp_maf, downcast="float", errors="coerce"
    )
    smaller_anno_df.gene_count = pd.to_numeric(
        smaller_anno_df.gene_count, downcast="integer", errors="coerce"
    )
    smaller_anno_df.dist_nearest_gene = pd.to_numeric(
        smaller_anno_df.dist_nearest_gene, downcast="float", errors="coerce"
    )
    smaller_anno_df.ld_buddies = pd.to_numeric(
        smaller_anno_df.ld_buddies, downcast="float", errors="coerce"
    )

    # remove +inf and -inf
    smaller_anno_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    keep_anno_df = smaller_anno_df[~smaller_anno_df.isna().any(1)].copy()
    del smaller_anno_df

    all_lead_snps = get_snps_to_match(input_snps_file)
    snps_to_match, snps_excluded = write_excluded_snps(
        all_lead_snps, keep_anno_df, OutObj.get("excluded_snps_file")
    )

    # calculate thresholds for max deviation from properties allowed
    thresholds = set_thresholds(
        max_dev_maf=0.05,
        max_dev_gene_count=0.50,
        max_dev_dist_nearest_gene=0.50,
        max_dev_friends_ld=0.50,
    )

    ###
    ###   start matching
    ###

    # start multithreading - one SNP per thread
    num_threads = cpu_count()
    mstart = time.time()
    logger.debug(
        "Generating {:,} control sets for each lead snp using {:,} cores.".format(
            n_matches, num_threads - 1
        )
    )

    pool = Pool(processes=4, maxtasksperchild=10)
    partial_get_matched_snps = partial(
        get_matched_snps,
        n_matches=n_matches,
        anno_df=keep_anno_df,
        thresholds=thresholds,
    )
    matched_snps_list = pool.map(partial_get_matched_snps, snps_to_match)

    pool.close()
    pool.join()
    logger.info(
        "Done matching SNPS. Took {:.2f} minutes.".format((time.time() - mstart) / 60)
    )

    # consolidate results
    snps_to_match_ordered, matched_snps, num_matched_found = zip(*matched_snps_list)
    column_names = ["Set_{}".format(x + 1) for x in range(n_matches)]
    agg_df = pd.DataFrame(
        np.vstack(matched_snps), columns=column_names
    )  # ADD 1 to start counting control sets at 1 instead of zero
    agg_df["snps_to_match"] = snps_to_match_ordered
    reorder_columns = ["snps_to_match"] + agg_df.columns[:-1].tolist()

    # remove input snps that didn't ahve any matched snps
    raw_final_matched_df = agg_df.loc[:, reorder_columns].copy()  # reorganize columns
    input_snp_wo_matched_snps = raw_final_matched_df[
        (raw_final_matched_df.loc[:, column_names] == "None").all(1)
    ].snps_to_match.tolist()

    ###
    ###   evaluate quality of matching
    ###

    match_qual_df = write_match_quality(
        snps_to_match_ordered,
        matched_snps,
        n_matches,
        num_matched_found,
        OutObj.get("quality_score_file"),
        fsummary,
    )
    fsummary.writelines("Done! Took {:.2f} minutes".format((time.time() - start) / 60))
    fsummary.close()

    # exclude snps with < 90% of the required number of control snps
    ex_insuff_input_snps = set(
        match_qual_df.loc[
            match_qual_df["n_matches_below_0.9"] == True, "lead_snp"
        ].values.tolist()
    )
    write_excluded_snps_w_insuff_matches(
        OutObj.get("excluded_snps_file"), ex_insuff_input_snps
    )

    ###
    ###    write final matched set
    ###
    snps_to_remove = set(input_snp_wo_matched_snps).union(ex_insuff_input_snps)
    logger.info(
        "* Excluded {:,} input snp(s) because of poor matching quality.".format(
            len(snps_to_remove)
        )
    )

    final_matched_df = raw_final_matched_df.loc[
        ~raw_final_matched_df["snps_to_match"].isin(snps_to_remove)
    ].copy()
    
    final_matched_df.rename(
        columns={"snps_to_match": "lead_snp"}, inplace=True
    )  # necessary for downstream analysis
    final_matched_df.to_csv(OutObj.get("matched_snps_file"), sep="\t", index=False)
    logger.debug("Wrote matched SNPs to: {}".format(OutObj.get("matched_snps_file")))

    return OutObj


# -----------
# MAIN
# -----------

if __name__ == "__main__":

    ###
    #   set up
    ###

    input_snps_file, n_matches, ld_buddies_r2, db_file, output_root = parse_input_args()

    _ = match_snps(input_snps_file, n_matches, ld_buddies_r2, db_file, output_root)
