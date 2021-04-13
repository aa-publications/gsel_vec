#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-03-06 10:39:18


import os
import pandas as pd
from datetime import datetime

import pickle

DATE = datetime.now().strftime("%Y-%m-%d")


###
###    PATHS # noqa
###

ANNO_DIR = "../data/anno_dict"  # noqa
OUTPUT_FILE = "../create_annotations"  # noqa
anno_path_dict = {
    "argweave": os.path.join(ANNO_DIR, "argweave_snpsnap_eur_ld0.1_collection.pickle"),
    "betascore": os.path.join(
        ANNO_DIR, "betascore_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "fst_eas_afr": os.path.join(
        ANNO_DIR, "fst_eas_afr_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "fst_eur_afr": os.path.join(
        ANNO_DIR, "fst_eur_afr_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "fst_eur_eas": os.path.join(
        ANNO_DIR, "fst_eur_eas_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "gerp": os.path.join(ANNO_DIR, "gerp_snpsnap_eur_ld0.1_collection.pickle"),
    "iES_Sabeti": os.path.join(
        ANNO_DIR, "iES_Sabeti_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "linsigh": os.path.join(ANNO_DIR, "linsigh_snpsnap_eur_ld0.1_collection.pickle"),
    "phastCon100": os.path.join(
        ANNO_DIR, "phastCon100_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "phyloP100": os.path.join(
        ANNO_DIR, "phyloP100_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "xpehh_afr2_eas": os.path.join(
        ANNO_DIR, "xpehh_afr2_eas_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "xpehh_afr2_eur": os.path.join(
        ANNO_DIR, "xpehh_afr2_eur_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "xpehh_eas_eur": os.path.join(
        ANNO_DIR, "xpehh_eas_eur_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "B2": os.path.join(
        ANNO_DIR, "B2_snpsnap_eur_ld0.1_collection.pickle"
    ),
    "geva_allele_age": os.path.join(
        ANNO_DIR, "geva_allele_age_snpsnap_eur_ld0.1_collection.pickle"
    ),
}

# -----------
# MAIN
# -----------
all_summary_df = pd.DataFrame()
for this_anno, anno_file in anno_path_dict.items():
    print(this_anno)
    anno = pickle.load(open(anno_file, "rb"))

    anno_df = pd.DataFrame.from_dict(anno, orient="index")
    anno_df.columns = ["value"]

    max = anno_df.max()[0]
    min = anno_df.min()[0]
    mean = anno_df.mean()[0]
    median = anno_df.median()[0]
    std = anno_df.std()[0]
    q25 = anno_df.quantile(0.25)[0]
    q75 = anno_df.quantile(0.75)[0]
    q95 = anno_df.quantile(0.95)[0]
    q5 = anno_df.quantile(0.05)[0]

    this_summary_df = pd.DataFrame(
        {
            "annotation": [this_anno],
            "max": [max],
            "min": [min],
            "mean": [mean],
            "median": [median],
            "std": [std],
            "q25": [q25],
            "q75": [q75],
            "q95": [q95],
            "q5": [q5],
        }
    )
    all_summary_df = all_summary_df.append(this_summary_df)


all_summary_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
