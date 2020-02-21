#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-02-16 12:55:00


import os
import sys
import gzip
import pickle
import numpy as np
import pandas as pd
from datetime import datetime

from time import time

DATE = datetime.now().strftime('%Y-%m-%d')


fstart=time()



import argparse

parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
parser.add_argument('annotation', action='store', type=str)
results = parser.parse_args()
annotation = results.annotation



# -----------
# PATH
# -----------

pop="eur"
dist="ld0.1"
snpsnap_db_file = f"/dors/capra_lab/projects/gwas_allele_age_evolution/data/snpsnap_database/raw_download/{pop}/{dist}_collection.tab.gz"


output_dir="/scratch/abraha1/gsel_/gsel_pipeline_vec/create_annotations/anno_dicts"

# annotaiton dictionary

param_dict={'argweave':{'anno_label':"argweave",
     'anno_suffix':".bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/argweaver/by_chr_tmrca_european",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/argweaver/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'betascore':{'anno_label':"beta_score",
     'anno_suffix':"_betas_v1.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/SNPsnap_BetaScores",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/SNPsnap_BetaScores/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},



     'fst_eas_afr':{'anno_label':"fst_eas_afr",
     'anno_suffix':"_fst_eas_afr.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eas_afr",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eas_afr/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'fst_eur_afr':{'anno_label':"fst_eur_afr",
     'anno_suffix':"_fst_eur_afr.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eur_afr",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eur_afr/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'fst_eur_eas':{'anno_label':"fst_eur_eas",
     'anno_suffix':"_fst_eur_eas.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eur_eas",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/fst/eur_eas/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'gerp':{'anno_label':"gerp",
     'anno_suffix':"_gerp_rs.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/GERP_RS/by_chr",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/GERP_RS/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'iES_Sabeti':{'anno_label':"iES_sabeti_eur",
     'anno_suffix':"_iES_Sabeti.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/iES_Sabeti_EUR",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/iES_Sabeti_EUR/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},



     'linsigh':{'anno_label':"linsight",
     'anno_suffix':"_LINSIGHT.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/linsight",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/linsight/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'phastCon100':{'anno_label':"phastcon100",
     'anno_suffix':".bed.gz" ,
     'annotation_dir':"/dors/capra_lab/data/evolutionary_conservation/phastcons/hg19_100way/by_chr",
     'snpsnap_anno_file':'/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/phastcon100/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz'},


     'phyloP100':{'anno_label':"phyloP100",
     'anno_suffix':".bed.gz" ,
     'annotation_dir':"/dors/capra_lab/data/evolutionary_conservation/phyloP/phyloP_100way_hg19/by_chr",
     'snpsnap_anno_file':'/dors/capra_lab/projects/gwas_allele_age_evolution/data/annotations/phyloP/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz'},



     'xpehh_afr2_eas':{'anno_label':"xpehh_afr2_eas",
     'anno_suffix':"_xpehh_afr2_eas.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/AFR2_EAS",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/AFR2_EAS/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'xpehh_afr2_eur':{'anno_label':"xpehh_afr2_eur",
     'anno_suffix':"_xpehh_afr2_eur.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/AFR2_EUR",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/AFR2_EUR/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},


     'xpehh_eas_eur':{'anno_label':"xpehh_eas_eur",
     'anno_suffix':"_xpehh_eas_eur.bed.gz" ,
     'annotation_dir':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/EAS_EUR",
     'snpsnap_anno_file':"/dors/capra_lab/projects/gwas_allele_age_evolution/data/iHS_ress/xpehh_results/xpehh_results_all/EAS_EUR/snpsnap_intersect/all_chr_snpsnap_intersect.bed.gz"},

     }


# -----------
# FUNCTIONS
# -----------
def check_file_paths(param_dict):

 missing_files = pd.DataFrame()
 for anno, anno_dict in param_dict.items():

     suffix = anno_dict['anno_suffix']
     annodir = anno_dict['annotation_dir']



     for num in np.arange(1,23,1):
         anno_path= os.path.join(annodir, "chr{}{}".format(num, suffix))
         if not os.path.exists(anno_path):
             temp_df = pd.DataFrame({'anno':[anno], 'path':[anno_path]})
             missing_files = missing_files.append(temp_df)


 return missing_files



# -----------
# MAIN
# -----------



# missing_file_df = check_file_paths(param_dict)
# print(f"number of missing files: {missing_file_df.shape[0]}")

# load snpsnap databset and record all the snps
sstart = time()
print(f"loading snpsnap db: {os.path.basename(snpsnap_db_file)}")
snpdb_df = pd.read_csv(snpsnap_db_file, sep="\t", usecols=['snpID'])
snpdb_snp = snpdb_df.snpID.unique()
print(f"Done loading {len(snpdb_snp):,} snpsnap snps. Took {(time()-sstart)/60:.2f} minutes.")



# loop through each annotation
# for this_anno in param_dict.keys():
this_anno = annotation


# annotation file
anno_time = time()
print(f">>loading {this_anno}...")
anno_file = param_dict[this_anno]['snpsnap_anno_file']
full_anno_df = pd.read_csv(anno_file, sep="\t", names=['chr','start','end','value'], usecols=[0,1,2,3])
print(f"Done loading {full_anno_df.shape[0]:,} snps for {this_anno}, took {(time()-anno_time)/60:.2f}.")


# full_anno_df = pd.DataFrame()
# for chromosome in np.arange(1,23,1):
#     print(chromosome)
#     anno_file=os.path.join(param_dict[this_anno]['annotation_dir'], "chr{}{}".format(chromosome, param_dict[this_anno]['anno_suffix']))
#     # load and format dict
#     anno_df = pd.read_csv(anno_file, sep="\t", names=['chr','start','end','value'])
#     full_anno_df = full_anno_df.append(anno_df)
#
# print(f"Done loading {full_anno_df.shape[0]:,} snps for {this_anno}, took {(time()-anno_time)/60:.2f}.")



full_anno_df['chr_label'], full_anno_df['chr_num'] = full_anno_df.chr.str.split("chr", 1).str
full_anno_df.start = full_anno_df.start.map(int)
full_anno_df['rsID']= full_anno_df.chr_num.map(str) + ":" + full_anno_df.start.map(str)


# make dict
keep_anno_df  = full_anno_df[full_anno_df.rsID.isin(snpdb_snp)].copy()
rsID_anno_dict = dict(zip(keep_anno_df.rsID, keep_anno_df.value))
print(f">>Made dictionary with {len(rsID_anno_dict):,}.")


# how many snps are commone between the two.. ?
dict_output_file = os.path.join(output_dir, f'{this_anno}_snpsnap_{pop}_{dist}_collection.pickle')
pickle.dump(rsID_anno_dict, open(dict_output_file, 'wb'))


print(f"Pickled dictionary. Total time: {(time()-anno_time)/60:.2f} minutes.")
