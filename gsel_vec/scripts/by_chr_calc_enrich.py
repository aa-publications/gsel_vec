#!/bin/python
# Abin Abraham
# created on: 2020-07-23 14:56:59
# calculate trait enrichment for each evolutionary annotation
#
#
#   * INPUTS:
#    - see argparse help
#    - if gsel_output_dir contains one folder per chromosome, set "mode" as "by_chr"; otherwise "all_chr"
#
#   * OUTPUT:
#    - tsv file with enrichments



import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import shutil
import time

import dask.dataframe as dd
from datetime import datetime
DATE = datetime.now().strftime("%Y-%m-%d")


# from .helper_dask_calc_enrich import *

## TODO: fix this path
# sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/run_mosaic_round_3_2021_01_27/2021_03_19_bolt-lmm_sstats")
# from helper_dask_calc_enrich import *
#


### ARGPARSE

if len(sys.argv) ==1:


    print("*** ")
    print("*** RUNNING IN DEBUG MODE ***")
    print("*** ")
    ANNO_FILE =Path("/Users/abinabraham/Documents/temp/repro_MENARCHE_AGE.sumstats")
    ANNO_SUMM_FILE=Path("/Users/abinabraham/Documents/temp/repro_MENARCHE_AGE.sumstats/chr1_repro_MENARCHE_AGE.sumstats/anno_genome_summary/genome_wide_summary_of_annotations.tsv")
    OUTPUT_DIR=Path("/Users/abinabraham/Documents/temp/repro_MENARCHE_AGE.sumstats")
    trait_name = 'test_trait'

else:
    import argparse
    parser = argparse.ArgumentParser(description="calc anno enrich")
    parser.add_argument("-g", dest="gsel_output_dir", action="store", type=str, help="full path to gsel_vec outputs for **ONE TRAIT**")

    results = parser.parse_args()
    ANNO_FILE = Path(results.gsel_output_dir.rstrip('\\'))
    ANNO_SUMM_FILE = list(Path(results.gsel_output_dir).glob('*/anno_genome_summary/genome_wide_summary_of_annotations.tsv'))[0] # pick one
    trait_name = Path(results.gsel_output_dir).stem # GSEL_Output directory is the name of trait_name
    OUTPUT_DIR=Path(results.gsel_output_dir)



print("---")
print(f"DATE: {DATE}")
print(f"script: {sys.argv[0]}")
print(f"gsel output dir: {ANNO_FILE}")
print(f"trait_name: {trait_name}")
print(f"OUTPUT_DIR: {OUTPUT_DIR}")
print("---")

# -----------
# FUNCTIONS
# -----------
def dask_load_annotation(trait_name, annotation, anno_files, cache_dir=None):

    dstart = time.time()
    dd_df = dd.read_csv(anno_files, sep="\t")
    n_lead_loci = dd_df.lead_snp.nunique().compute()

    # seperate into gwas and matched loci
    gwas_df = dd_df.loc[dd_df["set"] == "ld_snp", ["lead_snp", "anno"]].compute()
    matched_df = dd_df.loc[dd_df["set"] != "ld_snp", ["lead_snp", "set", "anno"]].compute()

    if cache_dir:
        ## cache
        gwas_out_file = cache_dir.joinpath(f"{trait_name}__{annotation}_raw_gwas_df.pickle")
        matched_out_file = cache_dir.joinpath(f"{trait_name}__{annotation}_raw_matched_df.pickle")
        gwas_df.to_pickle(gwas_out_file)
        matched_df.to_pickle(matched_out_file)



    print(f"Dask: Loaded {len(anno_files)} chromosomes in {took(time.time(), dstart )} minutes.")
    return gwas_df, matched_df

def calc_emp_pval(trait_summarized_df, matched_all_regions_df):

    # calc mean trait inside this function
    mean_trait = trait_summarized_df['anno'].mean()


    pstart = time.time()
    if mean_trait > np.median(matched_all_regions_df["anno"].values):
        emp_pval = np.sum(matched_all_regions_df["anno"].values >= mean_trait) / len(
            matched_all_regions_df["anno"].values
        )
    else:
        emp_pval = np.sum(matched_all_regions_df["anno"].values <= mean_trait) / len(
            matched_all_regions_df["anno"].values
        )
    # print(f"Done calc p-value. Took {took(time.time(), pstart)} mins.")
    return emp_pval

def count_before_after_na(matched_df,gwas_df):

    # remove na
    no_na_gwas_mean_df = gwas_df[~gwas_df.anno.isna()].copy()
    no_na_matched_df = matched_df[~matched_df.anno.isna()].copy()


    counts_after_summary_dict = {'n_og_matched_and_ld_snps': matched_df.shape[0],
                                 'n_no_na_matched_and_ld_snps': no_na_matched_df.shape[0],
                                 'n_matched_and_ld_snps_removed': matched_df.shape[0]-no_na_matched_df.shape[0],
                                 'percent_matched_and_ld_removed': ((matched_df.shape[0]-no_na_matched_df.shape[0])/matched_df.shape[0])*100,
                                 'n_og_trait_and_ld_snps': gwas_df.shape[0],
                                 'n_no_na_trait_and_ld_snps': no_na_gwas_mean_df.shape[0],
                                 'n_trait_and_ld_snps_removed': gwas_df.shape[0]-no_na_gwas_mean_df.shape[0],
                                 'percent_trait_and_ld_snps_removed': ((gwas_df.shape[0]-no_na_gwas_mean_df.shape[0])/gwas_df.shape[0])*100}


    return no_na_gwas_mean_df, no_na_matched_df, counts_after_summary_dict

def vec_sumarize_annotation(gwas_df, matched_df):


    no_na_gwas_mean_df, no_na_matched_df, counts_after_summary_dict = count_before_after_na(matched_df,gwas_df)
    # no_na_gwas_mean_df, no_na_matched_df, counts_after_summary_dict = count_before_after_na(subset_matched_df,subset_gwas_df)


    # trait associated regions --> take extreme value at each region
    no_na_gwas_mean_df['abs_anno'] = no_na_gwas_mean_df['anno'].abs()
    no_na_gwas_mean_df.sort_values(['lead_snp', 'abs_anno'], ascending=True, inplace=True)
    trait_summarized_df = no_na_gwas_mean_df[~no_na_gwas_mean_df.duplicated(subset=['lead_snp'], keep='last')].copy()
    del no_na_gwas_mean_df
    trait_summarized_df.drop(columns={'abs_anno'}, inplace=True)


    # matched regions --> take extreme value at each region
    no_na_matched_df['abs_anno'] = no_na_matched_df['anno'].abs()
    no_na_matched_df.sort_values(['lead_snp','set', 'abs_anno'], ascending=True, inplace=True)
    extreme_matched_mean_df = no_na_matched_df[~no_na_matched_df.duplicated(subset=['lead_snp','set'], keep='last')].copy()
    del no_na_matched_df
    extreme_matched_mean_df.drop(columns={'abs_anno'}, inplace=True)

    # --> take the mean for each matched set
    matched_all_regions_df = extreme_matched_mean_df.pivot(index='set', columns='lead_snp', values="anno").mean(1).reset_index()
    matched_all_regions_df.rename(columns={0:'anno'},inplace=True)


    # organize results
    descrip_stats_dict = descrip_stats_(trait_summarized_df, matched_all_regions_df)

    summarized_df_dict = {'trait_summarized_df':trait_summarized_df,
                        'extreme_matched_mean_df':extreme_matched_mean_df,
                        'matched_all_regions_df':matched_all_regions_df}




    return summarized_df_dict,counts_after_summary_dict, descrip_stats_dict

def descrip_stats_(trait_summarized_df, matched_all_regions_df):


    descrip_stats_dict = {}
    descrip_stats_dict['mean_trait_loci'] = trait_summarized_df['anno'].mean()
    descrip_stats_dict['median_trait_loci'] = trait_summarized_df.anno.quantile([0.5]).values[0]
    descrip_stats_dict['mean_matched_across_trait_loci'] = matched_all_regions_df.anno.values.mean()
    descrip_stats_dict['std_matched_across_trait_loci'] = matched_all_regions_df.anno.values.std()
    descrip_stats_dict['matched_5th']= matched_all_regions_df.anno.quantile(0.05)
    descrip_stats_dict['matched_25th']= matched_all_regions_df.anno.quantile(0.25)
    descrip_stats_dict['median_matched_loci']= matched_all_regions_df.anno.quantile(0.5)
    descrip_stats_dict['matched_75th']= matched_all_regions_df.anno.quantile(0.75)
    descrip_stats_dict['matched_95th']= matched_all_regions_df.anno.quantile(0.95)

    descrip_stats_dict['n_lead_loci_final'] = trait_summarized_df.lead_snp.nunique()
    descrip_stats_dict['n_matched_sets'] = matched_all_regions_df.set.nunique()

    return descrip_stats_dict

def calc_enrichment(genom_std_per_anno_dict, annotation, descrip_dict):


    genom_std = genom_std_per_anno_dict[annotation]
    enrich_per_mean_diff_by_genomstd = (descrip_dict['mean_trait_loci'] - descrip_dict['mean_matched_across_trait_loci']) / genom_std
    return enrich_per_mean_diff_by_genomstd

# %%

# -----------
# MAIN
# -----------
start = time.time()
lg = open(OUTPUT_DIR.joinpath(f'{trait_name}_enrich.log'), "w")
lg.write(f'{DATE}\n')
lg.write(f'{sys.argv[0]}\n')
lg.write(f"gsel output dir: {ANNO_FILE}\n")
lg.write(f"trait name: {trait_name}\n")
lg.write(f"output dir: {OUTPUT_DIR}\n")



# -----------
# load genome wide std of each evo annotation
# -----------
anno_genom_df = pd.read_csv(ANNO_SUMM_FILE, sep="\t")
genom_std_per_anno_dict = dict(zip(anno_genom_df["annotation"], anno_genom_df["std"]))

# %%
# -----------
# load annotation values
# -----------

# if by_chrm_bool:
# glob across (a) chromosomes and (b) selection measures
anno_file_list = list(ANNO_FILE.glob(f"*chr*/intermediate_analyses/selection_intersected_matched_sets/*/*_annotation_values.tsv"))


if len(anno_file_list) == 0:
    lg.write("Exiting. No annotation_values files to parse.")
    sys.exit("No annotation_values files to parse.")


lg.write(f"number of annotation files: {len(anno_file_list)}\n")


# %%
# -----------
# calc enrichment per annotation
# -----------
took = lambda x,y: f"{(x-y)/60:.2f}"



### cache directory to speed up subsequent analsyes
cache_dir = OUTPUT_DIR.joinpath(f'cached_{trait_name}')
cache_dir.mkdir(parents=True, exist_ok=True)

enrich_summary_df = pd.DataFrame()
uniq_annos = np.unique([anno_file.parts[-2] for anno_file in anno_file_list])
for n_anno, annotation in enumerate(uniq_annos):
    print(f"Anno {n_anno} out of {len(uniq_annos):,}")

    # select all chromosomes for this annotaiton
    anno_files = [anno_file for anno_file in anno_file_list if (anno_file.parts[-2] == annotation)]


    ### log
    anno_start = time.time()
    lg.write(f"\nParsing {annotation}. {n_anno} out of {len(uniq_annos):,}.\n")

    if len(anno_files) == 0:
        err = "No annotation_values files to parse\n"
        lg.write(err)
        sys.exit(err)
    elif len(anno_files) > 22:
        err = 'more than 22 files per annotation...\n'
        err = 'more than 22 files per annotation...\n'
        lg.write(err)
        sys.exit(err)
    else:
        lg.write(f"\tLoaded {len(anno_files)} chromosome annotation files.\n")


    ## load this annotation across all chromosomes
    gwas_df, matched_df = dask_load_annotation(trait_name, annotation, anno_files, cache_dir)

    # preflight check
    no_na_gwas_mean_df, no_na_matched_df, counts_after_summary_dict = count_before_after_na(matched_df,gwas_df)
    if (no_na_gwas_mean_df.shape[0] == 0) | (no_na_matched_df.shape[0]==0):

        id_dict = {'trait_name':trait_name,
                'annotation':annotation,
                "region_summary": 'extreme',
                "trait_summary": 'mean'}
        id_dict.update(counts_after_summary_dict)
        enrich_summary_df = enrich_summary_df.append(pd.DataFrame.from_records(id_dict, index=[0]))

    else:
        s = time.time()

        # summarize each locus (collapse over ld snps)
        summarized_df_dict, counts_after_summary_dict, descrip_dict = vec_sumarize_annotation(gwas_df, matched_df)

        # calc p-value
        emp_pval = calc_emp_pval(summarized_df_dict['trait_summarized_df'], summarized_df_dict['matched_all_regions_df'])

        # enrichment
        enrich_per_mean_diff_by_genomstd = calc_enrichment(genom_std_per_anno_dict, annotation, descrip_dict)

        if np.isnan(enrich_per_mean_diff_by_genomstd):
            emp_pvalue = np.nan

        enrich_results_dict = {'trait_name':trait_name,
                              'annotation':annotation,
                              'enrich_per_mean_diff_by_genomstd':enrich_per_mean_diff_by_genomstd,
                              'emp_pval':emp_pval,
                              'regions_summary':'extreme',
                              'trait_summary':'mean'}

        enrich_results_dict.update(descrip_dict)
        enrich_results_dict.update(counts_after_summary_dict)

        # save results
        enrich_summary_df = enrich_summary_df.append(pd.DataFrame.from_records(enrich_results_dict, index=[0]))

    print(f"Done with annotation. Took {took(time.time(), anno_start)} minutes.")

# if not cache_only:
    # write file with enrichments
summary_outfile = OUTPUT_DIR.joinpath(f"{trait_name}_extreme_regions_mean_enrichment_all_annotation.tsv")
enrich_summary_df.to_csv(summary_outfile, index=False, sep="\t", float_format='%.5f')
print(f"Done with enrichments. Took {took(time.time(), start )} minutes.")
lg.write(f"Done with enrichments. Took {took(time.time(), start )} minutes.")
lg.close()
