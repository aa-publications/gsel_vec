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

    # %%
    fig, axs = plt.subplots(nrows=len(anno_label_list),  sharex=True)
    for ind, annotation in enumerate(anno_label_list):

        ax = axs[ind]
        anno_df = anno_zscore_df.loc[anno_zscore_df['annotation']==annotation].copy()
        ax.plot(anno_df['xind'], anno_df['mean_controls'], 'o', color='gray', markersize=2)
        ax.plot(anno_df['xind'], anno_df['lead_snp_anno'], 'D', color='indianred', markersize=4, label=annotation.replace("_", " "))
        upper_std = anno_df['mean_controls']+anno_df['std_controls']
        lower_std = anno_df['mean_controls']-anno_df['std_controls']
        ax.plot([anno_df['xind'],anno_df['xind']], [lower_std, upper_std], '-', color='gray',)
        # ax.set_title(annotation)
        ax.set_xticks(np.arange(0, anno_df.shape[0]))
        ax.set_xticklabels([lead_xind_snp_dict[x] for x in np.arange(0, anno_df.shape[0])], rotation=270)
        # ax.legend(framealpha=0.5, borderpad=0, handletextpad=0.1, loc='upper right')
        ax.set_ylabel(annotation)

    # fig.text(0.5, 0.04, 'Lead SNPs', ha='center')
    # fig.text(0.04, 0.5, 'Annotation Value', va='center', rotation='vertical')

    # %%
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
    DIR="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/demov1/selection_intersected_matched_sets"
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
    # anno_pval_df, anno_zscore_df, _ = load_zscore_and_pvalue(intersectAnnoOutputObj, anno_path_dict )
    anno_pval_df, anno_zscore_df, _ = dev_load()


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

# %%
    cg = sns.clustermap(wide_zscore_df, row_cluster=False, cmap='vlag', center=0, annot=wide_pval_label_df, fmt='s')
    cg.ax_heatmap.set_yticklabels([x.get_text() for x in cg.ax_heatmap.yaxis.get_majorticklabels()], rotation=0)



    # %%
    return cg


    ###
    ###    plot actual annotation value with confidence intervals
    ###



def prep_radar_data(TraitEnrichOutObj):

    summary_enrich_file = TraitEnrichOutObj.get('trait_enrichment')
    summary_df = pd.read_csv(summary_enrich_file, sep="\t")
    plt_df = summary_df.copy()
    enrichs = summary_df['enrich_per_mean_diff_by_genomstd'].values.tolist()
    enrichs += enrichs[:1]


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

    return enrichs, angles, sig_angles, sig_enrichs, annotation_labels

def dev_prep_radar_data():

    # summary_enrich_file = TraitEnrichOutObj.get('trait_enrichment')
    # summary_enrich_file="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/demov1/trait_enrichments/trait_enrichment.tsv"
    summary_enrich_file="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/demov2/intermediate_analyses/trait_enrichments/trait_enrichment.tsv"
    summary_df = pd.read_csv(summary_enrich_file, sep="\t")
    plt_df = summary_df.copy()
    enrichs = summary_df['enrich_per_mean_diff_by_genomstd'].values.tolist()
    enrichs += enrichs[:1]


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

    return enrichs, angles, sig_angles, sig_enrichs, annotation_labels


def plot_radar( angles, annotation_labels, enrichs, sig_angles, sig_enrichs):


    # *****  DEV   *****
    # enrichs, angles, sig_angles, sig_enrichs, annotation_labels = dev_prep_radar_data()

    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(5,5))
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(annotation_labels, color='red', size=8,rotation=300)


    # set number of co-centric circles
    min_y = np.floor((np.min(enrichs)-0.1)*10)/10
    max_y = np.ceil((np.max(enrichs)+0.1)*10)/10
    max_ = np.max( [np.abs(max_y), np.abs(min_y)])


    ax.set_ylim(-max_, max_)
    ax.set_yticks(np.linspace(-max_, max_, 5))
    ax.set_yticklabels(["{:.2f}".format(x) for x in np.linspace(-max_, max_, 5)], ha='center')
    ax.set_rlabel_position(90,)


    ax.plot(angles, enrichs, linewidth=0, linestyle='solid', marker='o', markerfacecolor='gray', color='k', markersize=8)
    ax.plot(sig_angles, sig_enrichs, linewidth=0, linestyle='solid', marker='o',markerfacecolor='indianred',color='k', markersize=8, label='p<0.05')


    # relabel the radial axies
    lb_angles = np.linspace(0,2*np.pi,len(ax.get_xticklabels())+1)
    lb_angles[np.cos(lb_angles) < 0] = lb_angles[np.cos(lb_angles) < 0] + np.pi
    lb_angles = np.rad2deg(lb_angles)
    labels = []
    for label, angle in zip(ax.get_xticklabels(), lb_angles):
        x,y = label.get_position()
        lab = ax.text(x,y-0.17, label.get_text(), transform=label.get_transform(),
                      ha=label.get_ha(), va=label.get_va(), fontsize=10, color='black')
        lab.set_rotation(angle)
        labels.append(lab)

    ax.set_xticklabels([])
    # ax.annotate("enrichment",(-pi/18,0.1), fontsize=12, color='gray')

    # ax.annotate("enrichment:= [mean(gwas_loci) - mean(matched_loci)] / std(matched_loci) ",(0.35,0.01), fontsize=8, xycoords='figure fraction', color='indianred', annotation_clip=False, ha='center')

    gridlines = ax.yaxis.get_gridlines()
    # gridlines[ind].set_linewidth(2.5)
    # [x.set_color("black") for x in gridlines]
    ax.spines['polar'].set_visible(False)


    ax.legend( fontsize=8, loc='upper right')
    ax.fill(angles, enrichs, 'skyblue', alpha=0.3)
    ax.fill_between(np.linspace(0, 2*np.pi, 100), -0.001, 0.001, color='goldenrod', zorder=10) # <-- Added here
    [x.set_color('goldenrod') if x.get_position() == (0,0) else '' for x in ax.get_yticklabels()]

    # plt.subplots_adjust(top=0.8)
    # ax.annotate('Trait-Enrichment', xy=(0.5, 0.9), xycoords='figure fraction', fontsize='x-large', weight='bold', ha='center')
    # plt.savefig(test.pdf)
    # ax.set_title("Trait-Enrichment", loc='center')
    # plt.tight_layout()


# %%
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

    # make heatmap
    cg = plot_zscore_heatmap(intersectAnnoOutputObj, anno_path_dict)
    cg.savefig(OutObj.get('zscore_heatmap'))
    logger.info("Heatmap of trait-associated loci saved.")

    # make mean annotation value plot
    fig, axs = plot_mean_annotation_values(intersectAnnoOutputObj, anno_path_dict)
    fig.savefig(OutObj.get('mean_anno_scatter'))
    logger.info("Mean annotation value scatter plot saved.")


    # make radar plot
    enrichs, angles, sig_angles, sig_enrichs, annotation_labels = prep_radar_data(TraitEnrichOutObj)
    fig, ax = plot_radar(angles, annotation_labels, enrichs, sig_angles, sig_enrichs)
    fig.savefig(OutObj.get('trait_enrich_radar'))
# %%

