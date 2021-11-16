#!/usr/bin/env python



import os
import sys
import argparse
from pathlib import Path
from gsel_vec import run_evo_sig

ROOT = Path(os.path.abspath(os.path.dirname(__file__)))
TESTING_DATA = ROOT.parent.joinpath('gsel_vec', 'demo', 'input_data','bmi_small.test')



if __name__ =="__main__":

    parser = argparse.ArgumentParser(
        description="Get evolutionary signatures using GWAS summary stats."
    )
    parser.add_argument(
        "--analysis_name",
        action="store",
        type=str,
        help="name of this analysis no spaces allowed",
        default="testing",
    )
    parser.add_argument(
        "--gwas_summary_file",
        action="store",
        type=str,
        help="full path to gwas summary stats file",
        default=TESTING_DATA,
    )
    parser.add_argument(
        "--outputpath",
        action="store",
        type=str,
        help="full path to create output directory named 'analysis_name'",
        default=os.getcwd(),
    )

    args = parser.parse_args()
    analysis_name = args.analysis_name
    gwas_summary_file = args.gwas_summary_file
    outputpath = args.outputpath



    run_evo_sig.main(analysis_name, gwas_summary_file, outputpath)
    