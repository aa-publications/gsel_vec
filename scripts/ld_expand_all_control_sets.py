#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import argparse
import time
import math
import numpy as np
import pandas as pd
from datetime import datetime

import psutil
import gc
from functools import partial
from glob import glob

import logging

logger = logging.getLogger("main.{}".format(__name__))


from .helper_general import Outputs, report_mem

DATE = datetime.now().strftime("%Y-%m-%d")


# -----------
# FUNCTIONS
# -----------


def parse_input_args():

    # TODO: convert all args to requried args

    parser = argparse.ArgumentParser(description="Run expand control set.")

    parser.add_argument(
        "lead_ld_counts_file", action="store", type=str, help="table of lead GWAS SNPs"
    )

    parser.add_argument(
        "gwas_snps_r2_file",
        action="store",
        type=str,
        help="table of r2 between gwas snps",
    )

    parser.add_argument(
        "matched_file",
        action="store",
        type=str,
        help="input/lead gwas snps with matched SNPs",
    )

    parser.add_argument(
        "control_ld_dir",
        action="store",
        type=int,
        help="directory with ld SNPs for control snps",
    )

    parser.add_argument("output_dir", action="store", type=str, help="output_dir")

    # retrieve passed arguments
    args = parser.parse_args()
    lead_ld_counts_file = args.lead_ld_counts_file
    gwas_snps_r2_file = args.gwas_snps_r2_file
    matched_file = args.matched_file
    control_ld_dir = args.control_ld_dir
    output_dir = args.output_dir

    return (
        lead_ld_counts_file,
        gwas_snps_r2_file,
        matched_file,
        control_ld_dir,
        output_dir,
    )


def force_ld_control_snps(ld_df, controls_snps):
    """
    * does NOT ensure that all control snps are present in the ld _df
    """

    # keep all control snps in "SNP_A"
    temp_a = ld_df.loc[ld_df["SNP_A"].isin(controls_snps)].copy()
    # keep all control snps in "SNP_B"
    temp_b = ld_df.loc[ld_df["SNP_B"].isin(controls_snps)].copy()

    # TODO: remove excess column
    # flip SNP A and SNP B assignments
    temp_b.rename(columns={"SNP_A": "SNP_B", "SNP_B": "SNP_A"}, inplace=True)
    # temp_b.columns = ['CHR_B', 'BP_B', 'SNP_B', 'CHR_A', 'BP_A', 'SNP_A', 'R2', 'DP']
    temp_b = temp_b.loc[:, temp_a.columns]

    forced_ld_df = pd.concat([temp_a, temp_b])

    forced_ld_df.SNP_A = forced_ld_df.SNP_A.astype("category")
    forced_ld_df.SNP_B = forced_ld_df.SNP_B.astype("category")
    forced_ld_df.R2 = forced_ld_df.R2.astype("float16")

    logger.info(
        "* Note: Out of {:,} control snps, {:,} ({:.2f}%) snps were not found after running plink ld.".format(
            len(controls_snps),
            len(controls_snps) - forced_ld_df.SNP_A.nunique(),
            (len(controls_snps) - forced_ld_df.SNP_A.nunique()) / len(controls_snps),
        )
    )

    return forced_ld_df


def load_ld_df(LD_DIR):

    # TODO, low priority, remove glob pattern and use lsit of files`
    ld_df = pd.DataFrame()
    for ldfile in glob(LD_DIR + "/*.ld"):

        temp_ld_df = pd.read_csv(
            ldfile,
            sep="\s+",
            usecols=["SNP_A", "SNP_B", "R2"],
            dtype={"SNP_A": str, "SNP_B": np.str, "R2": np.float16},
        )
        temp_ld_df.SNP_A = temp_ld_df.SNP_A.astype("category")
        temp_ld_df.SNP_B = temp_ld_df.SNP_B.astype("category")

        ld_df = ld_df.append(temp_ld_df)

    ld_df.SNP_A = ld_df.SNP_A.astype("category")
    ld_df.SNP_B = ld_df.SNP_B.astype("category")

    return ld_df


def shuffle_and_list_r2(x, max_n):
    "return a list with all ld snps for a given control snp. NOTE: the control snp is appended as the first element."
    # note: max_n is the length of the list to be return;
    #       this will be equal to the dimension of the mask array storing the LD SNPs, which should include the control SNP

    num_ld_snps = max_n - 1  # since we will append the control snp

    control_snp = [x.SNP_A.values[0]]
    ld_snps_of_control_snp = x.SNP_B
    ldsnps_lst = list(ld_snps_of_control_snp)

    r2_of_ld_snps = x.R2
    r2_lst = list(r2_of_ld_snps)

    if len(ldsnps_lst) >= num_ld_snps:
        return [1] + r2_lst[
            0:num_ld_snps
        ]  # I checked to make sure the incides & slices make sense

    else:

        mulpt = int(
            np.ceil(num_ld_snps / len(ldsnps_lst))
        )  # I checked to make sure that number of times it is tiled is correct

        return [1] + (r2_lst * mulpt)[0:num_ld_snps]


def shuffle_and_list(x, max_n):
    "return a list with all ld snps for a given control snp. NOTE: the control snp is appended as the first element."
    # note: max_n is the length of the list to be return;
    #       this will be equal to the dimension of the mask array storing the LD SNPs, which should include the control SNP

    num_ld_snps = max_n - 1  # since we will append the control snp

    control_snp = [x.SNP_A.values[0]]
    ld_snps_of_control_snp = x.SNP_B
    ldsnps_lst = list(ld_snps_of_control_snp)

    # r2_of_ld_snps =  x.R2
    # r2_lst =  list(r2_of_ld_snps)

    if len(ldsnps_lst) >= num_ld_snps:
        return (
            control_snp + ldsnps_lst[0:num_ld_snps]
        )  # I checked to make sure the incides & slices make sense

    else:

        mulpt = int(
            np.ceil(num_ld_snps / len(ldsnps_lst))
        )  # I checked to make sure that number of times it is tiled is correct

        return control_snp + (ldsnps_lst * mulpt)[0:num_ld_snps]


def same_chrom_check(og_lead_control_ld_expanded_df):
    lead_control_ld_expanded_df = og_lead_control_ld_expanded_df.copy()
    lead_control_ld_expanded_df.fillna("None", inplace=True)
    ncontrol_set = lead_control_ld_expanded_df.shape[1] - 4
    temp_chrom_df = pd.concat(
        (
            lead_control_ld_expanded_df.loc[:, ["lead_snp", "R2", "lead_snp_bool"]],
            lead_control_ld_expanded_df.loc[
                :, [f"Set_{x}" for x in np.arange(1, ncontrol_set + 1)]
            ].applymap(lambda x: x.split(":")[0]),
        ),
        axis=1,
    )
    bool_df = (
        temp_chrom_df.loc[:, temp_chrom_df.columns.difference(["R2", "lead_snp_bool"])]
        .replace("None", np.nan)
        .groupby("lead_snp")
        .nunique()
        == 1
    )

    return bool_df, temp_chrom_df


def explode(df, lst_cols, fill_value="", preserve_index=False):

    # make sure `lst_cols` is list-alike
    if (
        lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))
    ):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = pd.DataFrame(
        {col: np.repeat(df[col].values, lens) for col in idx_cols}, index=idx
    ).assign(**{col: np.concatenate(df.loc[lens > 0, col].values) for col in lst_cols})
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = res.append(df.loc[lens == 0, idx_cols], sort=False).fillna(fill_value)
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res


def to_row_per_csnp(selected_ld_snps, n_control_snps):
    raw_control_df = pd.DataFrame(
        selected_ld_snps.swapaxes(1, 2).reshape(1, -1, n_control_snps).squeeze()
    )
    raw_control_df.replace("", np.nan, inplace=True)
    control_df = raw_control_df[~raw_control_df.isnull().any(1)].copy()
    control_df.reset_index(inplace=True, drop=True)
    control_df.columns = [
        "Set_{}".format(num) for num in np.arange(1, control_df.shape[1] + 1)
    ]
    return control_df


def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])


#
# def report_mem():
#     process = psutil.Process(os.getpid())
#
#     rss = convert_size(np.round(process.memory_info().rss))
#     uss = convert_size(np.round(process.memory_info().uss))
#
#     return f"{rss} rss,  {uss} uss"


def format_gwas_lead_and_ld_snps(gwas_lead_ld_df, matched_lead_snps):
    # ensure that every lead snp from 'matched_lead_snps'has at least one row where r2==1 with itself

    gwas_lead_ld_df["lead_snp_bool"] = (
        gwas_lead_ld_df["lead_snp"] == gwas_lead_ld_df["ld_snp"]
    )

    missing_lead_snps = list(
        set(matched_lead_snps)
        - set(
            gwas_lead_ld_df.loc[
                gwas_lead_ld_df["lead_snp_bool"] == True
            ].lead_snp.unique()
        )
    )
    temp_lead_df = pd.DataFrame(
        {
            "lead_snp": list(missing_lead_snps),
            "ld_snp": missing_lead_snps,
            "R2": 1.0,
            "lead_snp_bool": True,
        }
    )
    final_lead_ld_df = pd.concat((gwas_lead_ld_df, temp_lead_df), axis=0)

    final_lead_ld_df.sort_values(
        ["lead_snp", "lead_snp_bool", "R2"], ascending=False, inplace=True
    )
    final_lead_ld_df.reset_index(inplace=True, drop=True)
    final_lead_ld_df.R2.replace("NONE", np.nan, inplace=True)
    final_lead_ld_df.R2 = pd.to_numeric(final_lead_ld_df.R2)

    return final_lead_ld_df


def dev_load_control_ld_snps(chrom_num, control_snps):

    ld_structure_file = os.path.join(
        control_ld_dir, "ld_for_control_snps_chr{}.ld".format(chrom_num)
    )
    ld_struct_df = pd.read_csv(
        ld_structure_file,
        sep="\s+",
        usecols=["SNP_A", "SNP_B", "R2"],
        dtype={"SNP_A": str, "SNP_B": np.str, "R2": np.float16},
    )

    ld_struct_df.SNP_A = ld_struct_df.SNP_A.astype("category")
    ld_struct_df.SNP_B = ld_struct_df.SNP_B.astype("category")

    snps_in_chrom = set(
        [snp for snp in control_snps if (int(snp.split(":")[0]) == chrom_num)]
    )
    forced_ld_df = force_ld_control_snps(ld_struct_df, snps_in_chrom)
    forced_ld_df["ld_bin"] = 1 / np.round(
        forced_ld_df.R2 * 10
    )  # take the reciprocal so that it sorts in decreasing r2

    # SHUFFLE to simulate random selection later on...
    rand_ind = np.arange(
        forced_ld_df.shape[0]
    )  # I checked to make sure same # of rand_ind are generated as num rows
    np.random.seed(22)
    np.random.shuffle(rand_ind)
    forced_ld_df["rand_int"] = rand_ind
    forced_ld_df.sort_values(["SNP_A", "ld_bin", "rand_int"], inplace=True)
    forced_ld_df.drop(["rand_int"], axis=1, inplace=True)

    # consider these filters:
    # load ld pairs w/ *forced* the input/lead linked control_snp to be in teh first column
    # NOTE: * filter for r2 -->  we remove variants w/ r2 above the higher threhsold in the window
    #       * this means LD SNPs will be picked randomly from stratified LD bins that are decreasing in r2
    #       * ld_df will not have a row for lead SNPs w/o ld snps
    temp_filt_ld_df = forced_ld_df.loc[
        forced_ld_df["R2"] <= uppper_ld_threshold, ["SNP_A", "SNP_B", "R2"]
    ].copy()

    # remove aall SNP_A SNP_A r2=1 hits so that it doest not get tiled when there are not LD SNPs
    filt_r2_ld_df = temp_filt_ld_df.loc[
        ~temp_filt_ld_df.duplicated(["SNP_A", "SNP_B"], keep=False), :
    ].copy()
    filt_ld_df = temp_filt_ld_df.loc[
        ~temp_filt_ld_df.duplicated(["SNP_A", "SNP_B"], keep=False), ["SNP_A", "SNP_B"]
    ].copy()
    filt_ld_df.reset_index(drop=True, inplace=True)

    return forced_ld_df


def dev_filter_shuffle(forced_ld_df, control_snps, uppper_ld_threshold):

    # make sure ld_df already has control_snps in the first column ...

    # filter
    control_snp_mask = forced_ld_df.SNP_A.isin(control_snps)  # keep only control snps
    upper_r2_mask = (
        forced_ld_df.R2 <= uppper_ld_threshold
    )  # keep r2 values ≤ upper r2 threshold
    no_dup_snpa_mask = ~forced_ld_df.duplicated(
        ["SNP_A", "SNP_B"], keep=False
    )  # remove same snp in snp_A and snp_B with r2 ==1
    filt_ld_df = forced_ld_df.loc[
        control_snp_mask & upper_r2_mask & no_dup_snpa_mask
    ].copy()

    # create ld_bin
    filt_ld_df["ld_bin"] = 1 / np.round(
        filt_ld_df.R2 * 10
    )  # take the reciprocal so that it sorts in decreasing r2

    # shuffle within each control snp and ld_bin to simulate random selection later on
    rand_ind = np.arange(
        filt_ld_df.shape[0]
    )  # I checked to make sure same # of rand_ind are generated as num rows
    np.random.seed(22)
    np.random.shuffle(rand_ind)
    filt_ld_df["rand_int"] = rand_ind
    filt_ld_df.sort_values(["SNP_A", "ld_bin", "rand_int"], inplace=True)
    filt_ld_df.drop(["rand_int", "ld_bin"], axis=1, inplace=True)

    # add control snps that do not have any LD snps with it
    missing_control_snps = set(control_snps) - set(filt_ld_df.SNP_A.unique())
    miss_df = pd.DataFrame(
        {"SNP_A": list(missing_control_snps), "SNP_B": "NONE", "R2": "NONE"}
    )
    complete_ld_df = pd.concat((miss_df, filt_ld_df))

    complete_ld_df.SNP_A = complete_ld_df.SNP_A.astype("category")
    complete_ld_df.SNP_B = complete_ld_df.SNP_B.astype("category")

    return complete_ld_df, missing_control_snps


def array_to_per_control_snp_df(selected_ld_snps, n_control_snps):

    raw_control_df = pd.DataFrame(
        selected_ld_snps.transpose().reshape(-1, n_control_snps)
    )
    raw_control_df.drop(
        np.where(raw_control_df.isnull().all(1))[0], axis=0, inplace=True
    )
    raw_control_df.reset_index(inplace=True, drop=True)
    raw_control_df.columns = [
        "Set_{}".format(num) for num in np.arange(1, raw_control_df.shape[1] + 1)
    ]

    return raw_control_df


def ld_expand_all_control_snps(
    lead_ld_counts_file,
    gwas_snps_r2_file,
    matched_file,
    control_ld_dir,
    output_root,
    ld_thresholds=["ld<=1.0", "ld<=0.9"],
):
    """
    NOTE: ld_thresholds must be
                * a list with each element have this format: 'ld<={}'.format(r2_threshold); r2_threshold must be 1.0, 0.9, 0.8 etc..
                * the r2 threshold must go from high to low


    """

    tstart = time.time()

    ###
    ###    outputs
    ###

    logger.info("Starting to LD expand all control snps.")
    output_dir = os.path.join(output_root, "ld_expanded_control_and_input_snps")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj.add_output(
        "ld_expanded_output", "ld_expanded_all_control_sets.tsv", add_root=True
    )
    OutObj.add_output(
        "ld_r2_expanded_output", "r2_ld_expanded_all_control_sets.tsv", add_root=True
    )

    # temp file that will be out of order
    OutObj.add_output("temp_ld_exp_output", "temp_ld_exp.tsv", add_root=True)
    OutObj.add_output("temp_r2_ld_exp_output", "temp_r2_ld_exp.tsv", add_root=True)

    ###
    ###   load matched control sets
    ###

    matched_df = pd.read_csv(matched_file, sep="\t")
    control_snps = matched_df.iloc[:, 1:].values.flatten()
    n_lead_snps = matched_df.shape[0]
    n_control_snps = matched_df.shape[1] - 1
    matched_lead_snps = matched_df.lead_snp.unique()
    logger.info(
        "LD expanding control sets for {:,} lead snps with {:,} control snps.".format(
            n_lead_snps, n_control_snps
        )
    )

    ###
    ###   load lead and LD pairs of input snps
    ###

    # keep only snps with matched control snps (since some gwas lead snps may not be in the snpsnap database..)
    gwas_lead_ld_df = pd.read_csv(gwas_snps_r2_file, sep="\t")
    gwas_lead_ld_df.drop(
        np.where(~gwas_lead_ld_df.lead_snp.isin(matched_lead_snps))[0],
        axis=0,
        inplace=True,
    )

    if gwas_lead_ld_df.shape[0] == 0:
        sys.exit("No lead snps left to match after intersection with the snpsnap database.")


    # ensure every matched_lead_snp has at least a row with r2==1
    clean_lead_ld_df = format_gwas_lead_and_ld_snps(gwas_lead_ld_df, matched_lead_snps)
    del gwas_lead_ld_df

    # define the order of lead snps
    ordered_lead_snps = clean_lead_ld_df[
        ~clean_lead_ld_df.duplicated(subset=["lead_snp"], keep="first")
    ].lead_snp.tolist()
    ordered_lead_snps_dict = dict(
        zip(ordered_lead_snps, np.arange(len(ordered_lead_snps)))
    )

    ###
    ###    load the # of LD snps required for each lead gwas snp
    ###

    # read in number of LD snps required per lead gwas snps
    all_ld_counts_df = pd.read_csv(lead_ld_counts_file, sep="\t")

    # ensure that only the matched_lead_snps are present in the ld counts table
    all_ld_counts_df.drop(
        all_ld_counts_df[~all_ld_counts_df.lead_snp.isin(matched_df.lead_snp)].index,
        axis=0,
        inplace=True,
    )

    # reorder to match
    all_ld_counts_df.set_index("lead_snp", inplace=True)
    ordered_ld_counts_df = all_ld_counts_df.reindex(ordered_lead_snps)
    del all_ld_counts_df

    ###
    ###   load ld pairs for the **control snps** and *shuffle*
    ###

    ld_df = load_ld_df(control_ld_dir)
    ld_df.SNP_A = ld_df.SNP_A.astype("category")
    ld_df.SNP_B = ld_df.SNP_B.astype("category")

    # force the control snp to be in the first column
    forced_ld_df = force_ld_control_snps(
        ld_df, control_snps
    )  # note: not all lead snps are guarenteed to be present in forced_ld_df, AND that is FINE!
    del ld_df
    gc.collect()

    ###
    ###    loop through each ld bin
    ###

    # ld_bins_expanded_df = pd.DataFrame()
    # ld_bins_r2_expanded_df = pd.DataFrame()
    for ld_threshold in ld_thresholds:

        uppper_ld_threshold = float(ld_threshold.split("=")[1])
        lower_ld_threshold = uppper_ld_threshold - 0.1
        logger.info(
            "\tLD expanding: {} < r2 ≤ {}".format(lower_ld_threshold, uppper_ld_threshold)
        )
        logger.debug(f"Memory used before LD expanding: {report_mem()}")

        # -----------
        #  LD mask  detail:
        # -----------

        #  create a mask that will select SNPs in LD with the control snps
        #       * inner most array will be a vector of 1 or 0, serving as the mask to select LD SNPs for a specific control snp
        #       * (the first element will always be 1 to include the control SNP)
        #       * in numpy shape order ->
        #               * # of input/lead snps,
        #               * # of control snps requested,
        #               * length == max # of LD snps required out of all input snps
        #
        #
        #  Example:
        #           # control snp 1        # control snp 2
        #   [   [   [1, 1, 0, 0, 1, 0],   [1, 1, 1, 1, 0, 0], ...   ],         --> input/lead gwas snp #1
        #       [   [1, 1, 0, 0, 1, 0],   [1, 1, 1, 1, 0, 0], ...   ], ...  ]  --> input/lead gwas snp #2

        # -----------
        # create mask
        # -----------

        # create a binary vector for each input/lead snp where 1 is a placeholder to a required LD SNP
        # every lead snps will have at least one '1', which stands for the matched control snp
        # inner_vec_length = ordered_ld_counts_df.loc[:, ld_threshold].max()
        # padded_mask_per_input_snp = np.array( [  np.pad(np.ones(x), (0, inner_vec_length - (x)) , mode='constant') for x in ordered_ld_counts_df.loc[:, ld_threshold]])
        # assert padded_mask_per_input_snp.shape[0] == n_lead_snps, "number of lead snps don't match with mask array "

        # OLD VERSION WITH +1
        inner_vec_length = (
            ordered_ld_counts_df.loc[:, ld_threshold].max() + 1
        )  # +1 so that the control snp is included as the first element
        padded_mask_per_input_snp = np.array(
            [
                np.pad(np.ones(x + 1), (0, inner_vec_length - (x + 1)), mode="constant")
                for x in ordered_ld_counts_df.loc[:, ld_threshold]
            ]
        )
        assert (
            padded_mask_per_input_snp.shape[0] == n_lead_snps
        ), "number of lead snps don't match with mask array "

        # tile equal to number of control sets and swap axes to ensure input snps get the outer axes
        ld_mask = np.swapaxes(
            np.tile(padded_mask_per_input_snp, (n_control_snps, 1, 1)), 0, 1
        )
        assert (
            ld_mask.shape[1] == n_control_snps and ld_mask.shape[0] == n_lead_snps
        ), "array shapes don't match up with inputs"

        # -----------
        # create nested ld array
        # -----------

        # follow the same structure as teh nested_mask array
        # except inner most array will have LD SNPs instead of 1 or 0
        # note that each inner array should tile (this is to simulate sampling with replacement)
        # the length of the inner vector will be determined by the length of the innerrmost array of nested_mask

        # load ld pairs w/ *forced* the input/lead linked control_snp to be in teh first column
        # NOTE: * filter for r2 -->  we remove variants w/ r2 above the higher threhsold in the window
        #       * this means LD SNPs will be picked randomly from stratified LD bins that are decreasing in r2
        #       * ld_df will not have a row for lead SNPs w/o ld snps

        par_shuf = partial(shuffle_and_list, max_n=ld_mask.shape[2])
        par_r2_shuf = partial(shuffle_and_list_r2, max_n=ld_mask.shape[2])

        # for each input/lead snp
        for ind_snp, input_snp in enumerate(ordered_lead_snps):
            # set trace ...

            logger.debug(
                f"LD expandeing control snps for {ind_snp} out of {len(ordered_lead_snps)} lead snps."
            ) if ((ind_snp % 15) == 0) else None

            # get all control snps and remember the order!
            csnps_for_isnp = matched_df.loc[
                matched_df["lead_snp"] == input_snp,
                ["Set_{}".format(x) for x in range(1, n_control_snps + 1)],
            ].values.flatten()
            uniq_csnps = np.unique(csnps_for_isnp)
            csnps_indices = dict(
                zip(
                    uniq_csnps,
                    [
                        np.where(csnps_for_isnp == uniq_csnp)[0].tolist()
                        for uniq_csnp in uniq_csnps
                    ],
                )
            )

            # filter and group ld snps for each control snp - every control snp has a row
            # control snps w/o LD snps will have None
            complete_ld_df, miss_control_snps = dev_filter_shuffle(
                forced_ld_df, uniq_csnps, uppper_ld_threshold
            )
            ld_array_df = (
                complete_ld_df.groupby("SNP_A", sort=False)
                .apply(par_shuf)
                .reset_index()
            )
            ld_r2_array_df = (
                complete_ld_df.groupby("SNP_A", sort=False)
                .apply(par_r2_shuf)
                .reset_index()
            )

            ld_array_df.rename({0: "ld_array"}, axis=1, inplace=True)
            ld_r2_array_df.rename({0: "r2_ld_array"}, axis=1, inplace=True)

            # order control snps, this is necessary to preserve order between different ld bin iterations
            ld_array_df["order"] = ld_array_df.SNP_A.map(csnps_indices)
            ordered_ld_array_df = explode(
                ld_array_df, "order"
            )  # need to do this since a control snp could be present multiple times
            ordered_ld_array_df.sort_values("order", inplace=True)

            # order for the r2 dataframe as well
            ld_r2_array_df["order"] = ld_r2_array_df.SNP_A.map(csnps_indices)
            ordered_ld_r2_array_df = explode(ld_r2_array_df, "order")
            ordered_ld_r2_array_df.sort_values("order", inplace=True)

            # extract arrays of ld snps for each control snpd
            ldsnps_for_each_control_snp = ordered_ld_array_df.loc[
                :, "ld_array"
            ].values.tolist()
            r2_ldsnps_for_each_control_snp = ordered_ld_r2_array_df.loc[
                :, "r2_ld_array"
            ].values.tolist()

            # select LD SNPs based on ld_mask
            ld_snps_array = np.array(
                ldsnps_for_each_control_snp, dtype="object", order="F"
            )
            ld_snps_r2_array = np.array(
                r2_ldsnps_for_each_control_snp, dtype="object", order="F"
            )

            selected_ld_snps = np.where(
                np.array(ld_mask[ind_snp], dtype=bool), ld_snps_array, np.nan
            )
            selected_ld_r2_snps = np.where(
                np.array(ld_mask[ind_snp], dtype=bool), ld_snps_r2_array, np.nan
            )

            # convert to one column per control snp, one row per SNP in LD with control snp
            raw_control_df = array_to_per_control_snp_df(
                selected_ld_snps, n_control_snps
            )
            raw_r2_control_df = array_to_per_control_snp_df(
                selected_ld_r2_snps, n_control_snps
            )

            del ordered_ld_array_df, ordered_ld_r2_array_df, ld_r2_array_df, ld_array_df, selected_ld_snps, selected_ld_r2_snps
            gc.collect()
            # logger.debug(f"On {input_snp} ({ind_snp}) after append , using {report_mem()}")
            logger.debug(
                f"\t{report_mem()}: after randomly picking LD snps for control snps for this input snp. "
            )

            #
            # get the input snp and its ld snps column ready
            #

            input_snp_bool = clean_lead_ld_df["lead_snp"] == input_snp
            r2_threshold_bool = (clean_lead_ld_df["R2"] > lower_ld_threshold) & (
                clean_lead_ld_df["R2"] <= uppper_ld_threshold
            )

            input_lead_ld_df = clean_lead_ld_df.loc[
                (input_snp_bool & r2_threshold_bool)
            ].reset_index(drop=True)

            if uppper_ld_threshold != 1:
                temp_df = pd.DataFrame(
                    {
                        "lead_snp": [input_snp],
                        "ld_snp": [input_snp],
                        "R2": [1.0],
                        "lead_snp_bool": [True],
                    }
                )
                input_lead_ld_df = input_lead_ld_df.append(temp_df)

            # if the input snp, doesn't have any LD SNPS at this threshold, then ADD the r2 of the snp with itself
            # this can occur in the first bin (0.9-1.0), if criteria is not met
            if input_lead_ld_df.shape[0] == 0:
                input_lead_ld_df = pd.DataFrame(
                    {
                        "lead_snp": [input_snp],
                        "ld_snp": [input_snp],
                        "R2": [1.0],
                        "lead_snp_bool": [True],
                    }
                )

            assert (
                raw_control_df.shape[0] == input_lead_ld_df.shape[0]
            ), "Number of rows do not match between the two dataframes (lead snps and control snps) that are going to be concatentated "
            assert (
                raw_r2_control_df.shape[0] == input_lead_ld_df.shape[0]
            ), "Number of rows do not match between the two dataframes (lead snps and control snps) that are going to be concatentated "

            input_lead_ld_df.reset_index(drop=True, inplace=True)
            raw_control_df.reset_index(drop=True, inplace=True)
            raw_r2_control_df.reset_index(drop=True, inplace=True)
            lead_control_ld_expanded_df = pd.concat(
                (input_lead_ld_df, raw_control_df), axis=1, ignore_index=False
            )
            lead_control_ld_expanded_r2_df = pd.concat(
                (input_lead_ld_df, raw_r2_control_df), axis=1, ignore_index=False
            )
            del input_lead_ld_df, raw_control_df, raw_r2_control_df

            logger.debug(
                f"\t{report_mem()}: before writing control and r2 values for this input SNP."
            )
            # write

            temp_ld_file = OutObj.get("temp_ld_exp_output")
            with open(temp_ld_file, "a") as f:
                lead_control_ld_expanded_df.to_csv(
                    f, index=False, sep="\t", header=f.tell() == 0
                )

            temp_r2_ld_file = OutObj.get("temp_r2_ld_exp_output")
            with open(temp_r2_ld_file, "a") as f:
                lead_control_ld_expanded_r2_df.to_csv(
                    f, index=False, sep="\t", header=f.tell() == 0
                )

            del lead_control_ld_expanded_r2_df, lead_control_ld_expanded_df
            gc.collect()

    logger.debug(f"Done with all LD bins. Using {report_mem()}")

    full_df = pd.read_csv(OutObj.get("temp_ld_exp_output"), sep="\t")
    full_r2_df = pd.read_csv(OutObj.get("temp_r2_ld_exp_output"), sep="\t")

    logger.debug(f"Done reading in all ld expanded data. Using {report_mem()}")

    dedup_df = full_df[
        ~full_df.duplicated(
            subset=["lead_snp", "ld_snp", "R2", "lead_snp_bool"], keep="first"
        )
    ].sort_values(["lead_snp", "lead_snp_bool", "R2"], ascending=False)
    dedup_r2_df = full_r2_df[
        ~full_r2_df.duplicated(
            subset=["lead_snp", "ld_snp", "R2", "lead_snp_bool"], keep="first"
        )
    ].sort_values(["lead_snp", "lead_snp_bool", "R2"], ascending=False)

    dedup_df.replace("NONE", np.nan, inplace=True)
    dedup_r2_df.replace("NONE", np.nan, inplace=True)
    logger.debug(f"Done replacing NONE with nans. Using {report_mem()}")

    # write
    dedup_df.to_csv(
        OutObj.get("ld_expanded_output"), sep="\t", index=False, na_rep="None"
    )
    dedup_r2_df.to_csv(
        OutObj.get("ld_r2_expanded_output"), sep="\t", index=False, na_rep="None"
    )
    logger.info(
        "Wrote expanded control snps to {}".format(OutObj.get("ld_expanded_output"))
    )
    logger.info(
        "Done LD expanding control snps in {:.2f} minutes".format(
            (time.time() - tstart) / 60
        )
    )

    return OutObj


# -----------
# MAIN
# -----------


if __name__ == "__main__":

    lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, output_root = (
        parse_input_args()
    )
    ld_expand_all_control_snps(
        lead_ld_counts_file,
        gwas_snps_r2_file,
        matched_file,
        control_ld_dir,
        output_root,
    )
