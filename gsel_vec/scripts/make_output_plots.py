#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-01-07 14:04:55


import os
import sys
import numpy as np
import pandas as pd
from math import pi
from datetime import datetime
import time
from .helper_general import Outputs

import logging
logger = logging.getLogger("main.{}".format(__name__))


DATE = datetime.now().strftime('%Y-%m-%d')

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# %%
# -----------
# FUNCTIONS
# -----------

def set_up_outputs(OutputObj):

    OutputObj.add_output(
        "zscore_heatmap",
        "{}_zscore_heatmap.pdf".format(DATE), add_root=True)
    OutputObj.add_output(
        "mean_anno_scatter",
        "{}_mean_anno_scatter.pdf".format(DATE), add_root=True)
    OutputObj.add_output(
        "trait_enrich_radar",
        "{}_trait_enrich_radar.pdf".format(DATE), add_root=True)

    return OutputObj


def label_pvalue(x):

    if x >= 0.05:
        return ''
    elif (x < 0.05):
        return '*'
    else:
        return np.nan


def add_xindex_per_lead_snp(raw_anno_zscore_df):

    anno_zscore_df = raw_anno_zscore_df.copy()
    xind_anno_zscore_df  = anno_zscore_df.copy()
    xind_anno_zscore_df.drop_duplicates(subset=['lead_snp'], inplace=True)

    xind_anno_zscore_df['chr'] = xind_anno_zscore_df['lead_snp'].apply(lambda x: x.split(':')[0])
    xind_anno_zscore_df['pos'] = xind_anno_zscore_df['lead_snp'].apply(lambda x: x.split(':')[1])
    xind_anno_zscore_df.sort_values(['chr','pos'],ascending=True, inplace=True)

    xind_anno_zscore_df['xind'] = np.arange(0,xind_anno_zscore_df.shape[0])
    lead_snp_xind_dict = dict(zip(xind_anno_zscore_df['lead_snp'], xind_anno_zscore_df['xind']))
    lead_xind_snp_dict = dict(zip(xind_anno_zscore_df['xind'], xind_anno_zscore_df['lead_snp']))

    anno_zscore_df['xind'] = anno_zscore_df['lead_snp'].map(lead_snp_xind_dict)

    return anno_zscore_df, lead_xind_snp_dict


def plot_mean_annotation_values(intersectAnnoOutputObj, anno_path_dict):


    _, anno_zscore_df, anno_label_list = load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict)
    # _, anno_zscore_df, anno_label_list =  dev_load()
    anno_zscore_df, lead_xind_snp_dict = add_xindex_per_lead_snp(anno_zscore_df)


    fig, axs =plt.subplots(nrows=anno_zscore_df['annotation'].nunique(), sharex=True)

    annotations = anno_zscore_df['annotation'].unique()
    for ind, (ax,anno) in enumerate(zip(axs, annotations)):

        this_df = anno_zscore_df.loc[anno_zscore_df['annotation']==anno]
        ax.scatter(this_df['lead_snp'], this_df['lead_snp_anno'], marker='D', color='indianred', s=50, label='Mean')
        ax.errorbar(this_df['lead_snp'], this_df['mean_controls'], yerr=this_df['std_controls'], linewidth=1, fmt ='k+', label='matched set\n(mean+/-s.d.)')
        ax.set_ylabel(anno)

    axs[-1].set_xlabel('GWAS Loci')
    axs[0].legend(bbox_to_anchor=(1, 0.8))

    axs[0].set_title("Evolutionary Annotation Mean Value")
    _ = plt.xticks(rotation='vertical')

    return fig, axs


def load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict):
    anno_label_list, anno_paths = zip(*anno_path_dict.items())
    anno_pval_df = pd.DataFrame()
    anno_zscore_df = pd.DataFrame()
    for anno_label in anno_label_list:
        pval_file_label = "{}_pvalue".format(anno_label)
        zscore_file_label = "{}_zscore".format(anno_label)


        pval_file = intersectAnnoOutputObj.get(pval_file_label)
        zscore_file = intersectAnnoOutputObj.get(zscore_file_label)


        pval_df = pd.read_csv(pval_file, sep="\t")
        zscore_df = pd.read_csv(zscore_file, sep="\t")

        anno_label = anno_label.replace("_", " ")
        pval_df['annotation'] = anno_label
        zscore_df['annotation'] = anno_label

        anno_pval_df = anno_pval_df.append(pval_df)
        anno_zscore_df = anno_zscore_df.append(zscore_df)
    return anno_pval_df, anno_zscore_df, anno_label_list



def dev_load():

    from glob import glob
    # DIR="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/demov1/selection_intersected_matched_sets"
    DIR="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/bmi_small/intermediate_analyses/selection_intersected_matched_sets"
    annotaitons_dirs =  glob(DIR+"/*")


    anno_pval_df = pd.DataFrame()
    anno_zscore_df = pd.DataFrame()
    anno_label_list = []
    for anno_file in annotaitons_dirs:


        anno_label = os.path.basename(anno_file)
        pval_file_label = "{}_pvalue".format(anno_label)
        zscore_file_label = "{}_zscore".format(anno_label)


        pval_file = os.path.join(DIR,anno_label,  f"pvalue_{anno_label}_per_lead_snp.tsv")
        zscore_file = os.path.join(DIR, anno_label, f"zscore_{anno_label}_per_lead_snp.tsv")


        pval_df = pd.read_csv(pval_file, sep="\t")
        zscore_df = pd.read_csv(zscore_file, sep="\t")

        anno_label = anno_label.replace("_", " ")
        pval_df['annotation'] = anno_label
        zscore_df['annotation'] = anno_label

        anno_pval_df = anno_pval_df.append(pval_df)
        anno_zscore_df = anno_zscore_df.append(zscore_df)
        anno_label_list.append(anno_label)


    return anno_pval_df, anno_zscore_df, anno_label_list



def plot_zscore_heatmap(intersectAnnoOutputObj, anno_path_dict):


    ###
    ### load z-score and p-value for each annotation
    ###
    anno_pval_df, anno_zscore_df, _ = load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict )
    # import ipdb; ipdb.set_trace()
    # anno_pval_df, anno_zscore_df, _ = dev_load()


    # ### DEV
    # anno_pval_df, anno_zscore_df, _  = dev_load()


    ###
    ###    format for heatmap plot
    ###

    anno_pval_df['pvalue_label'] = ''
    anno_pval_df.loc[anno_pval_df['corrected_pval_benj_hoch'] <0.05, 'pvalue_label'] = '***'
    anno_pval_df.loc[anno_pval_df['pvalue'] <0.05, 'pvalue_label'] = '*'

    wide_pval_label_df = anno_pval_df.loc[:, ['lead_snp','pvalue_label','annotation']].pivot(index='annotation', columns='lead_snp', values='pvalue_label')
    wide_zscore_df = anno_zscore_df.loc[:, ['lead_snp','z_score','annotation']].pivot(index='annotation', columns='lead_snp', values='z_score')

    # convert nans to 0 to enable clustering
    wide_zscore_df[wide_zscore_df.isna()] = 0


    cg = sns.clustermap(wide_zscore_df, row_cluster=False, cmap='vlag', center=0, annot=wide_pval_label_df, fmt='s')
    cg.ax_heatmap.set_yticklabels([x.get_text() for x in cg.ax_heatmap.yaxis.get_majorticklabels()], rotation=0)



    return cg


    ###
    ###    plot actual annotation value with confidence intervals
    ###



def prep_radar_data(TraitEnrichOutObj, summary_enrich_file=False):

    if summary_enrich_file:
        summary_df = pd.read_csv(summary_enrich_file, sep="\t")
        summary_df.rename(columns={'n_lead_loci_final':'n_lead_loci', 'emp_pval':'emp_pvalue'}, inplace=True)
    else:
        summary_enrich_file = TraitEnrichOutObj.get('trait_enrichment')
        summary_df = pd.read_csv(summary_enrich_file, sep="\t")


    plt_df = summary_df.copy()
    enrichs = summary_df['enrich_per_mean_diff_by_genomstd'].values.tolist()
    enrichs += enrichs[:1]
    mean_n_lead_loci = np.round(summary_df['n_lead_loci'].mean())

    pvals = []
    for ind, row in plt_df.iterrows():
        this_pvalue = row.emp_pvalue < 0.05
        pvals.append(this_pvalue)

    pvals += pvals[:1]

    n_cat = summary_df.annotation.nunique()
    angles = [n / float(n_cat) * 2 * pi for n in range(n_cat)]
    angles += angles[:1]

    sig_angles = [angles[i]  for i,x in enumerate(pvals) if x]
    sig_enrichs = [enrichs[i]  for i,x in enumerate(pvals) if x]

    annotation_labels = summary_df.annotation.unique().tolist()
    anno_label_dict={'argweave': 'ARGweaver\n(TMRCA)',
                     'geva_allele_age': 'GEVA\n(allele age)',
                    'betascore': 'Beta Score',
                    'B2': 'Beta Score 2',
                    'linsigh': 'LINSIGHT',
                    'phastCon100': 'PhastCons',
                    'phyloP100': 'PhyloP',
                    'gerp': 'GERP',
                    'fst_eas_afr': 'Fst\neas-afr',
                    'fst_eur_afr': 'Fst\neur-afr',
                    'fst_eur_eas': 'Fst\neur-eas',
                    'iES_Sabeti': 'iES(Sabeti)',
                    'xpehh_afr2_eas': 'XP-EHH\nafr-eas',
                    'xpehh_afr2_eur': 'XP-EHH\nafr-eur',
                    'xpehh_eas_eur': 'XP-EHH\neas-eur'}

    annotation_labels = [anno_label_dict[x] for x in annotation_labels]

    return enrichs, angles, sig_angles, sig_enrichs, annotation_labels, mean_n_lead_loci

def plot_radar( angles, annotation_labels, enrichs, sig_angles, sig_enrichs, mean_n_lead_loci):


    # *****  DEV   *****
    # enrichs, angles, sig_angles, sig_enrichs, annotation_labels, mean_n_lead_loci = dev_prep_radar_data()
    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(10,10))
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(annotation_labels, color='dimgray', size=16,rotation=300)


    # set number of co-centric circles
    min_y = np.floor((np.min(enrichs)-0.1)*10)/10
    max_y = np.ceil((np.max(enrichs)+0.1)*10)/10
    max_ = np.max( [np.abs(max_y), np.abs(min_y)])


    ax.set_ylim(-max_, max_)
    ax.set_yticks(np.linspace(-max_, max_, 5))
    ax.set_yticklabels(["{:.2f}".format(x) for x in np.linspace(-max_, max_, 5)], ha='center',color='dimgray', size=16)
    ax.set_rlabel_position(90,)


    ax.plot(angles, enrichs, linewidth=0, linestyle='solid', marker='o', markerfacecolor='gray', color='k', markersize=12)
    ax.plot(sig_angles, sig_enrichs, linewidth=0, linestyle='solid', marker='o',markerfacecolor='indianred',color='k', markersize=12, label='p<0.05')

    gridlines = ax.yaxis.get_gridlines()
    # gridlines[ind].set_linewidth(2.5)
    # [x.set_color("black") for x in gridlines]
    ax.spines['polar'].set_visible(False)


    ax.legend( fontsize=16, loc='upper right')
    ax.fill(angles, enrichs, 'skyblue', alpha=0.3)
    ax.fill_between(np.linspace(0, 2*np.pi, 100), -0.001, 0.001, color='goldenrod', zorder=10) # <-- Added here
    [x.set_color('goldenrod') if x.get_position() == (0,0) else '' for x in ax.get_yticklabels()]


    return fig, ax


def make_output_plots(intersectAnnoOutputObj, TraitEnrichOutObj, anno_path_dict, output_root):


    start = time.time()
    logger.info("Generating final plots.")

    ###
    ### set up outputs
    ###

    output_dir = os.path.join(output_root, "plot_outputs")
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj = set_up_outputs(OutObj)

    if intersectAnnoOutputObj:
        # make heatmap
        cg = plot_zscore_heatmap(intersectAnnoOutputObj, anno_path_dict)
        cg.savefig(OutObj.get('zscore_heatmap'))
        logger.debug("[status] Heatmap of trait-associated loci saved.")
        plt.clf()

        # make mean annotation value plot
        # fig, axs = plot_mean_annotation_values(intersectAnnoOutputObj, anno_path_dict)
        # plt.tight_layout()
        # fig.savefig(OutObj.get('mean_anno_scatter'))
        # logger.debug("[status] Mean annotation value scatter plot saved.")
        # plt.clf()

    # make radar plot
    if TraitEnrichOutObj:
        enrichs, angles, sig_angles, sig_enrichs, annotation_labels, mean_n_lead_loci = prep_radar_data(TraitEnrichOutObj)
        fig, ax = plot_radar(angles, annotation_labels, enrichs, sig_angles, sig_enrichs, mean_n_lead_loci)
        fig.savefig(OutObj.get('trait_enrich_radar'))
        logger.debug("[status] Radar plot saved.")


# %%

