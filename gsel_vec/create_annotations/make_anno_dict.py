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
from pathlib import Path
from time import time

DATE = datetime.now().strftime('%Y-%m-%d')


OUTPUT_DIR = os.getcwd()




import argparse

parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
parser.add_argument('annotation_file_path', action='store', type=str)
results = parser.parse_args()
annotation_file_path = Path(results.annotation_file_path)

annotation_file_path.name

# -----------
# PATH
# -----------
anno_df = pd.read_csv( annotation_file_path, sep="\t")
anno_dict  = dict(zip(anno_df['chr_pos'], anno_df['value']))


output_pickle_file = os.path.join(OUTPTU_DIR, "{annotation_file_path.name}.pickle")
pickle.dump(output_pickle_file open(anno_dict, 'wb'))
print("pickled dictionary created ")

