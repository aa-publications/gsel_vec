#!/usr/bin/env python



import os
import sys
import argparse
from pathlib import Path
from gsel_vec import run_evo_sig
import pkg_resources


ROOT = Path(os.path.abspath(os.path.dirname(__file__)))
# TESTING_DATA = ROOT.parent.joinpath('gsel_vec', 'demo', 'input_data','bmi_small.test')

TESTING_DIR = pkg_resources.resource_filename("gsel_vec", "demo")
TESTING_DATA= os.path.join(TESTING_DIR, 'input_data','bmi_small.test')
print("data path si {}".format(TESTING_DATA))




if __name__ =="__main__":

    parser = argparse.ArgumentParser(
        description="Get evolutionary signatures using GWAS summary stats."
    )
    parser.add_argument(
        "-a", "--analysis_name",
        action="store",
        type=str,
        help="name of this analysis no spaces allowed",
        default="testing",
    )
    parser.add_argument(
        "-g", "--gwas_summary_file",
        action="store",
        type=str,
        help="full path to gwas summary stats file",
        default=TESTING_DATA,
    )
    parser.add_argument(
        "-o", "--output_path",
        action="store",
        type=str,
        help="full path to create output directory named 'analysis_name'",
        default=os.getcwd(),
    )

    args = parser.parse_args()
    analysis_name = args.analysis_name
    gwas_summary_file = args.gwas_summary_file
    outputpath = args.output_path



    run_evo_sig.main(analysis_name, gwas_summary_file, outputpath, ROOT)
