#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-18 09:38:36


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')


# PATH 
INPUT_FILE="test/input_data/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt"
root, basename = os.path.split(INPUT_FILE)
OUTPUT_FILE= os.path.join(root, 'formatted_'+basename)


df = pd.read_csv(INPUT_FILE, sep="\s+")


df.rename({'CHR':'chr',
           'POS':'pos',
           'SNP':'rsID',
           'P':'p'}, axis=1, inplace=True)


og_shape = df.shape

df = df[~df.chr.isnull()].copy()
df = df[~df.pos.isnull()].copy()
df = df[~df.p.isnull()].copy()
print("After removing null values, {:,} rows out of {:,} remain".format(df.shape[0], og_shape[0]))


df.chr = df.chr.astype(int)
df.pos = df.pos.astype(int)

new_df = df.loc[:, ['chr','pos','rsID','p']].copy()
new_df.to_csv(OUTPUT_FILE, sep="\t", index=False)