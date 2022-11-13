#!/usr/bin/env python
#
#
#
#
#
#
# Abin Abraham
# created on: 2020-01-05 08:20:56



import logging
import os
import sys
import time
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd
import pkg_resources
from pathlib import Path

from gsel_vec.scripts.check_ld_expanded_control_sets import check_ld_expanded_sets
from gsel_vec.scripts.clump_snps import clump_snp_list, clump_snps
from gsel_vec.scripts.combine_control_sets import combine_control_sets
from gsel_vec.scripts.get_ldsnps_for_control_snps import get_ldsnps_for_control_snps
from gsel_vec.scripts.helper_general import safe_mkdir, start_logger
from gsel_vec.scripts.intersect_annotation import intersect_all_annotations
from gsel_vec.scripts.ld_expand_all_control_sets import ld_expand_all_control_snps
from gsel_vec.scripts.match_snps import match_snps
from gsel_vec.scripts.helper_calc_genome_distribution_of_annotations import (
    calc_genome_distribution_of_annotations,
)
from gsel_vec.scripts.calc_trait_enrichment import calc_trait_entrichment
from gsel_vec.scripts.organize_final_outputs import organize_final_outputs
from gsel_vec.scripts.make_output_plots import make_output_plots

master_start = time.time()


# -----------
# Inputs and Outputs
# -----------
def main(analysis_name, gwas_summary_file, outputpath, ROOT, run_by_chr=False):


    outputdir = os.path.join(outputpath, analysis_name)
    anno_summary_dir = os.path.join(outputdir, "anno_genome_summary")
    print(f"Note: Will overwrite existing output files.")
    print(f"Outputs saved to: {outputdir}")

    # -----------
    # PIPELINE PARAMTERS
    # -----------

    # num_control_sets:             number of MAF and LD matched control snps
    # lead_snp_min_gwas_pvalue:     SNPs with GWAS p-value lower than this will be considered to identify potential lead snps
    # ld_snps_min_gwas_pvalue:      only SNPs with GWAS p-value lower than this will be considered to be in LD with the lead SNPs
    # min_r2_to_clump:              the r2 threshold used to identify independent GWAS regions (snps with r2 greater than this will be clumped together)
    # min_kb_from_index_to_clump:   the kb threshold used to identify independent GWAS regions
    # ld_expand_lead_snp_min_r2:    for matching, this is the r2 threshold to which we ld expand the lead snp
    # ldbuds_r2_threshold:          for snpsnap, the r2 threshold used to determine the number of LD buddies; must be formatted as 'friends_ld09' or 'friends_ld08' etc.
    # control_snps_ld_expand_r2:    r2 threshold to ld expand the control snps
    # ld_thresholds:                ld threshold to create matched control regions (must be in the form "ld<=1.0" or ["ld<=1.0","ld<=0.9"] if expanding to 0.8<r2≤1.
    # summary_type:                 how to summarize the evolutionary annotation over genomic regions (accepts min, max, median, mean)

    if analysis_name =="testing":
        num_control_sets = 50
    else:
        num_control_sets = 5000

    lead_snp_min_gwas_pvalue = 0.00000005
    ld_snps_min_gwas_pvalue = 0.00000005
    min_r2_to_clump = 0.1
    min_kb_from_index_to_clump = 500
    ld_expand_lead_snp_min_r2 = 0.9

    ldbuds_r2_threshold = "friends_ld09"
    control_snps_ld_expand_r2 = 0.9
    ld_thresholds = ["ld<=1.0"]
    summary_type = "max"


    # -----------
    # DEPENDENCIES
    # -----------


    # data paths
    ROOT=Path(ROOT)
    data_path = ROOT.parent.joinpath( 'data')


    snpsnap_db_file = os.path.join(data_path, "snpsnap_database/snpsnap_database_ld0.9.tab.gz")
    thous_gen_file = os.path.join(data_path, "1kg/EUR.chr{}.phase3.nodups")
    anno_dir = os.path.join(data_path, "anno_dict")
    anno_descrip_file = os.path.join(data_path, "anno_dict/anno_description.txt")

    descrip_df = pd.read_csv(anno_descrip_file, sep=",", names=['anno','filename','two_sided_bool'])
    descrip_df['path'] = descrip_df['filename'].apply(lambda x: os.path.join(anno_dir, x))
    anno_path_dict  = dict(zip(descrip_df['anno'], descrip_df['path']))
    two_sided_bool_dict = dict(zip(descrip_df['anno'], descrip_df['two_sided_bool']))


    # anno_path_dict = {
        # "argweave": os.path.join(anno_dir, "argweave_snpsnap_eur_ld0.1_collection.pickle"),
        # "betascore": os.path.join(
        #     anno_dir, "betascore_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "fst_eas_afr": os.path.join(
        #     anno_dir, "fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "fst_eur_afr": os.path.join(
        #     anno_dir, "fst_eur_afr_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "fst_eur_eas": os.path.join(
        #     anno_dir, "fst_eur_eas_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "gerp": os.path.join(anno_dir, "gerp_snpsnap_eur_ld0.1_collection.pickle"),
        # "iES_Sabeti": os.path.join(
        #     anno_dir, "iES_Sabeti_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "linsigh": os.path.join(anno_dir, "linsigh_snpsnap_eur_ld0.1_collection.pickle"),
        # "phastCon100": os.path.join(
        #     anno_dir, "phastCon100_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "phyloP100": os.path.join(
        #     anno_dir, "phyloP100_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "xpehh_afr2_eas": os.path.join(
        #     anno_dir, "xpehh_afr2_eas_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "xpehh_afr2_eur": os.path.join(
        #     anno_dir, "xpehh_afr2_eur_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "xpehh_eas_eur": os.path.join(
        #     anno_dir, "xpehh_eas_eur_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "geva_allele_age": os.path.join(
            # anno_dir, "geva_allele_age_snpsnap_eur_ld0.1_collection.pickle"
        # ),
        # "B2": os.path.join(anno_dir, "B2_snpsnap_eur_ld0.1_collection.pickle"),
    # }
    #
    # two_sided_bool_dict = {
    #     "argweave": True,
        # "betascore": False,
        # "fst_eas_afr": False,
        # "fst_eur_afr": False,
        # "fst_eur_eas": False,
        # "gerp": False,
        # "iES_Sabeti": False,
        # "linsigh": False,
        # "phastCon100": False,
        # "phyloP100": True,
        # "xpehh_afr2_eas": False,
        # "xpehh_afr2_eur": False,
        # "xpehh_eas_eur": False,
        # "xpehh_eas_eur": False,
        # "geva_allele_age": True,
        # "B2": False,
    }

#todo
# add cleaner name for plot

    # -----------
    # START LOGGER
    # -----------

    # TO DO: exit program if output dir already exists. ask for a different name
    safe_mkdir(outputdir)

    intermediate_dir = os.path.join(outputdir, "intermediate_analyses")
    final_output_dir = os.path.join(outputdir, "final_outputs")

    safe_mkdir(intermediate_dir)
    safe_mkdir(final_output_dir)

    logfile = os.path.join(outputdir, f"{analysis_name}.log")
    logger = start_logger(logfile)
    logger.debug(
        f"""

    Settings for this run:
    ---- Inputs/Outputs ---
    gwas summary file: {gwas_summary_file}
    output directory: {outputpath}

    ---- Parameters ---
    num_control_sets: {num_control_sets}
    lead_snp_min_gwas_pvalue: {lead_snp_min_gwas_pvalue}
    ld_snps_min_gwas_pvalue: {ld_snps_min_gwas_pvalue}
    min_r2_to_clump: {min_r2_to_clump}
    min_kb_from_index_to_clump: {min_kb_from_index_to_clump}
    ld_expand_lead_snp_min_r2: {ld_expand_lead_snp_min_r2}
    ldbuds_r2_threshold: {ldbuds_r2_threshold}
    control_snps_ld_expand_r2: {control_snps_ld_expand_r2}
    ld_thresholds: {ld_thresholds}
    summary_type: {summary_type}

    """
    )


    # -----------
    # START PIPELINE
    # -----------

    # clump snps based on LD
    OutObj = clump_snps(
        gwas_summary_file,
        intermediate_dir,
        thous_gen_file,
        lead_snp_min_gwas_pvalue=lead_snp_min_gwas_pvalue,
        ld_snps_min_gwas_pvalue=ld_snps_min_gwas_pvalue,
        min_r2_to_clump=min_r2_to_clump,
        min_kb_from_index_to_clump=min_kb_from_index_to_clump,
        ld_expand_lead_snp_min_r2=ld_expand_lead_snp_min_r2,
    )


    # match snps
    match_OutObj = match_snps(
        OutObj.get("lead_snp_ld_pairs_r2"),
        num_control_sets,
        ldbuds_r2_threshold,
        snpsnap_db_file,
        intermediate_dir,
    )


    # get ld snps for control snps
    csnps_file = match_OutObj.get("matched_snps_file")
    ldexp_OutObj = get_ldsnps_for_control_snps(
        csnps_file, thous_gen_file, intermediate_dir, control_snps_ld_expand_r2
    )


    # ld expand control snps for each input/lead snp
    lead_snps_ld_counts_file = OutObj.get("bin_by_ld_file")
    gwas_snps_r2_file = OutObj.get("lead_snp_ld_pairs_r2")
    matched_file = match_OutObj.get("matched_snps_file")
    control_ld_dir = ldexp_OutObj.get("ld_snps_for_control_snps_dir")

    ldexp_match_OutObj = ld_expand_all_control_snps(
        lead_snps_ld_counts_file,
        gwas_snps_r2_file,
        matched_file,
        control_ld_dir,
        intermediate_dir,
        ld_thresholds=ld_thresholds,
    )


    # check ld expanded control sets
    ld_expanded_control_sets_file = ldexp_match_OutObj.get("ld_expanded_output")
    ld_expanded_control_sets_r2_file = ldexp_match_OutObj.get("ld_r2_expanded_output")
    (
        match_summary_by_params_df,
        ldscore_lead_and_ld_df,
        match_quality_per_lead_snp_df,
    ) = check_ld_expanded_sets(
        snpsnap_db_file,
        ld_expanded_control_sets_file,
        lead_snps_ld_counts_file,
        ld_expanded_control_sets_r2_file,
        ldbuds_r2_threshold,
        intermediate_dir,
    )


    # intersect annotation
    intersectAnnoOutputObj = intersect_all_annotations(
        anno_path_dict,
        two_sided_bool_dict,
        summary_type,
        ld_expanded_control_sets_file,
        intermediate_dir,
    )

    # calculate genome wide summary of annotations
    anno_genom_summary_file = calc_genome_distribution_of_annotations(
        anno_path_dict, anno_summary_dir
    )

    if run_by_chr:
        print("skipping trait-wide enrichments")
        # make heatmap of one chromosome
        intersect_ouputs = make_output_plots(
            intersectAnnoOutputObj, None, anno_path_dict, final_output_dir
        )

    else:

        # calculate trait-wide enrichment
        TraitEnrichOutObj = calc_trait_entrichment(
            intersectAnnoOutputObj, anno_genom_summary_file, anno_path_dict, final_output_dir
        )


        # make final plots
        intersect_ouputs = make_output_plots(
            intersectAnnoOutputObj, TraitEnrichOutObj, anno_path_dict, final_output_dir
        )




    # organize final outputs
    finalOutObj = organize_final_outputs(
        intersectAnnoOutputObj,
        anno_path_dict,
        match_quality_per_lead_snp_df,
        match_summary_by_params_df,
        ldbuds_r2_threshold,
        final_output_dir,
    )



    logger.debug(
        "[status] Done with master.py. Took {:.2f} minutes.".format(
            (time.time() - master_start) / 60
        )
    )


# %%
if __name__ == "__main__":
    main()
