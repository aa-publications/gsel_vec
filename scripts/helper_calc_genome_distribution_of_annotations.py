#!/bin/python
# This script will ...
#
#
#
#       TO DO: annotation paths are hard coded and if two-sided or one-sided p-values are also hard coded
#
#
# Abin Abraham
# created on: 2020-01-05 08:20:56

import os
import sys
import numpy as np
import pandas as pd
import pickle 
from .helper_general import safe_mkdir



from datetime import datetime
DATE = datetime.now().strftime('%Y-%m-%d')


import logging
logger = logging.getLogger("main.{}".format(__name__))


# -----------
# FUNCTIONS
# -----------

# %%
# -----------
# MAIN
# -----------

def calc_genome_distribution_of_annotations(anno_path_dict, anno_summary_dir):

    logger.info(f"Calculating the genome-wide summary of evolutionary annotations.")


    safe_mkdir(anno_summary_dir)


    all_summary_df = pd.DataFrame()
    for this_anno, anno_file in anno_path_dict.items():

        logger.debug("Calculating genome summary for {}.".format(this_anno))
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

    anno_genom_summary_file = os.path.join(anno_summary_dir, 'genome_wide_summary_of_annotations.tsv')
    all_summary_df.to_csv(anno_genom_summary_file, sep="\t", index=False)
    logger.info(f"Wrote the genome-wide summary of evolutionary annotations.")

    return anno_genom_summary_file
