#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import pickle
import time

import numpy as np
import pandas as pd

from scipy import stats
from .helper_general import Outputs

import logging
logger = logging.getLogger("main.{}".format(__name__))

from datetime import datetime
DATE = datetime.now().strftime("%Y-%m-%d")



# -----------
# FUNCTIONS
# -----------

def set_up_outputs(OutObj):
    # within that dir,  set up ouput files
    OutObj.add_output(
        "trait_enrichment",
        "trait_enrichment.tsv",
        add_root=True)

    return OutObj


def get_extreme(arr):

    return arr[np.argmax(np.abs(arr))]



def load_annotation(anno_file):

    # load one annotation
    anno_df = pd.read_csv(anno_file, sep="\t")
    n_lead_loci = anno_df.lead_snp.nunique()

    # seperate into gwas and matched loci
    gwas_df = anno_df.loc[anno_df["set"] == "ld_snp", ["lead_snp", "anno"]].copy()
    matched_df = anno_df.loc[
        anno_df["set"] != "ld_snp", ["lead_snp", "set", "anno"]
    ].copy()

    return gwas_df, matched_df



def summarize_annotation(
    gwas_df, matched_df, region_summary=get_extreme, trait_summary=np.mean
):

    # remove Na's
    no_na_gwas_mean_df = gwas_df[~gwas_df.anno.isna()].copy()
    no_na_matched_df = matched_df[~matched_df.anno.isna()].copy()


    logger.debug(
        "Removed {:,}({:.2f} trait-associated snps.".format(gwas_df.shape[0]-no_na_gwas_mean_df.shape[0],(gwas_df.shape[0]-no_na_gwas_mean_df.shape[0])/gwas_df.shape[0]))

    logger.debug( "Removed {:,}({:.2f}%) control-associated snps.".format(matched_df.shape[0]-no_na_matched_df.shape[0], (matched_df.shape[0]-no_na_matched_df.shape[0])/matched_df.shape[0]))


    # summarize selection annotation
    # trait-associated regions: summarize lead and LD SNPs
    logger.debug(
        "{:,} unique trait-associated regions identified.".format(no_na_gwas_mean_df.lead_snp.nunique()))

    trait_summarized_df = (
        no_na_gwas_mean_df.groupby("lead_snp")["anno"]
        .apply(lambda x: region_summary(x.values))
        .reset_index()
    )

    # matched regions: summarize lead and LD SNPs
    logger.debug(
        "{:,} matched regions identified per trait-associated regions.".format(no_na_matched_df.set.nunique()))

    matched_mean_df = (
        no_na_matched_df.groupby(["lead_snp", "set"])["anno"]
        .apply(lambda x: region_summary(x.values))
        .reset_index()
    )
    # summarize across all trait-associated regions
    matched_all_regions_df = (
        matched_mean_df.groupby("set")
        .apply(lambda x: trait_summary(x))["anno"]
        .reset_index()
    )

    return trait_summarized_df, matched_all_regions_df, matched_mean_df

def calc_emp_pval(mean_trait, matched_all_regions_df):
    if mean_trait > np.median(matched_all_regions_df["anno"].values):
        emp_pval = np.sum(matched_all_regions_df["anno"].values >= mean_trait) / len(
            matched_all_regions_df["anno"].values
        )
    else:
        emp_pval = np.sum(matched_all_regions_df["anno"].values <= mean_trait) / len(
            matched_all_regions_df["anno"].values
        )

    return emp_pval

def append_enrich_row_for_anno(annotation, emp_pval, region_summary, trait_summary, trait_summarized_df, matched_all_regions_df, summary_df):
    mean_trait_loci = trait_summarized_df.anno.values.mean()
    median_trait_loci = trait_summarized_df.anno.quantile([0.5]).values[0]
    mean_matched_across_trait_loci = matched_all_regions_df.anno.values.mean()
    std_matched_across_trait_loci = matched_all_regions_df.anno.values.std()
    matched_5th, median_matched_loci, matched_95th = matched_all_regions_df.anno.quantile(
        [0.05, 0.5, 0.95]
    )

    summary_df = summary_df.append(
        pd.DataFrame(
            {
                "annotation": [annotation],
                "mean_trait_loci": [mean_trait_loci],
                "median_trait_loci": [median_trait_loci],
                "emp_pvalue": [emp_pval],
                "mean_matched_across_trait_loci": [mean_matched_across_trait_loci],
                "std_matched_across_trait_loci": [std_matched_across_trait_loci],
                "matched_5th": [matched_5th],
                "median_matched_loci": [median_matched_loci],
                "matched_95th": [matched_95th],
                "region_summary": [region_summary[0]],
                "trait_summary": [trait_summary[0]],
                "n_lead_loci": [trait_summarized_df.lead_snp.nunique()],
                "n_matched_loci": [matched_all_regions_df.set.nunique()],
            }
        )
    )
    return summary_df

def calc_trait_entrichment(intersectAnnoOutputObj, anno_genom_summary_file, anno_path_dict, output_root):

    t_start = time.time()
    logger.info("Calculating trait enrichment.")


    ###
    ### load genome wide annotation summary for each evolutionary measure
    ###

    anno_label_list, anno_paths = zip(*anno_path_dict.items())
    anno_genom_summary_df = pd.read_csv(anno_genom_summary_file, sep="\t")

    ###
    ### set up outputs
    ###

    output_dir = os.path.join(output_root, "trait_enrichments")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)


    ###
    ### how to summarize at trait and regions
    ###

    trait_summary = ("mean", np.mean)
    region_summary = ("extreme", get_extreme)


    ###
    ### for each annotation, calculate trait-wide summaries
    ###

    summary_df = pd.DataFrame()
    for anno_label in anno_label_list:

        anno_file = intersectAnnoOutputObj.get("{}_annotation_values".format(anno_label))
        if not os.path.isfile(anno_file):
            logger.info("Trait enrichment was not calculate because could not find annotation file for {}.".format(anno_label))
            continue

        # load and summarize_annotation across all trait-assoicated loci
        gwas_df, matched_df = load_annotation(anno_file)
        trait_summarized_df, matched_all_regions_df, matched_mean_df = summarize_annotation(
            gwas_df,
            matched_df,
            region_summary=region_summary[1],
            trait_summary=trait_summary[1],
        )

        # calc one-sided p-value
        mean_trait = trait_summarized_df.mean().values[0]
        emp_pval = calc_emp_pval(mean_trait, matched_all_regions_df)
        summary_df = append_enrich_row_for_anno(anno_label, emp_pval, region_summary, trait_summary, trait_summarized_df, matched_all_regions_df, summary_df)





    # add genome_wide summary column
    summary_df["genome_wide_anno_std"] = summary_df["annotation"].map(
        dict(zip(anno_genom_summary_df["annotation"], anno_genom_summary_df["std"]))
    )


    ###
    ### calculate trait-wide enrichments
    ###

    summary_df["enrich_per_mean_diff_by_genomstd"] = (
        summary_df["mean_trait_loci"] - summary_df["mean_matched_across_trait_loci"]
    ) / summary_df["genome_wide_anno_std"]




    summary_df.to_csv(OutObj.get('trait_enrichment'), index=False, sep="\t")
    logger.info("Done with trait enrichment. Took {:.2f} minutes.".format( (time.time()- t_start)/60))
    logger.debug("[status] Done; outputs written to: {}".format(OutObj.get('trait_enrichment')))

    return OutObj