#!/usr/bin/env python



import os
import sys
import argparse
from pathlib import Path
from gsel_vec import run_evo_sig
from gsel_vec.scripts.by_chr_calc_enrich import func_by_chr_calc_enrich
import pkg_resources


# ROOT = Path(os.path.abspath(os.path.dirname(__file__)))
# print("ROOT path is {}".format(ROOT))

TESTING_DIR = pkg_resources.resource_filename("gsel_vec", "demo")
TESTING_DATA= os.path.join(TESTING_DIR, 'input_data','bmi_small.test')

DATA_ROOT_DIR = pkg_resources.resource_filename("gsel_vec", "data")
print("data path is {}".format(DATA_ROOT_DIR))


if __name__ =="__main__":


    # -----------
    # argparse
    # -----------

    ## if no input arguments -> run testing schema

    if len(sys.argv) ==1:


        print("*** ")
        print("*** RUNNING IN DEBUG MODE ***")
        print("*** ")

        analysis_name = 'testing'
        gwas_summary_file = TESTING_DATA
        outputpath = os.getcwd()
        run_by_chr=False
        run_by_chr_enrich = False

    else:
        parser = argparse.ArgumentParser(
            description="Get evolutionary signatures using GWAS summary stats."
        )
        parser.add_argument(
            "-a", "--analysis_name",
            action="store",
            type=str,
            help="name of this analysis no spaces allowed",
            required=True,
        )
        parser.add_argument(
            "-g", "--gwas_summary_file",
            action="store",
            type=str,
            help="full path to gwas summary stats file",
            default=None,
        )
        parser.add_argument(
            "-o", "--output_path",
            action="store",
            type=str,
            help="full path to create output directory named 'analysis_name'",
            required=True,
        )

        parser.add_argument(
        '-r', '--run_by_chr',
            action='store_true',
            default=False,
            help='flag whether this is just one chromosome',
        )

        parser.add_argument(
        '-e', '--run_by_chr_enrich',
            action='store_true',
            default=False,
            help='flag whether to run trait enrichment after runnign all chromosomes seperately',
        )


        args = parser.parse_args()
        analysis_name = args.analysis_name
        gwas_summary_file = args.gwas_summary_file
        outputpath = args.output_path
        run_by_chr_enrich = args.run_by_chr_enrich
        if args.run_by_chr:
            run_by_chr=True
        else:
            run_by_chr=False

    # -----------
    # run pipeline
    # -----------

    if run_by_chr_enrich:

        if run_by_chr:
            print("flags are not properly specified. Pick either run_by_chr and run_by_chr_enrich flag")

        else:
            outputdir = os.path.join(outputpath, analysis_name)
            _ = func_by_chr_calc_enrich(outputdir)

    # default
    # this will run either the entire pipeline or just one chromosome.
    else:
        run_evo_sig.main(analysis_name, gwas_summary_file, outputpath, DATA_ROOT_DIR, run_by_chr)
