#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
import logging
import psutil
import math

DATE = datetime.now().strftime("%Y-%m-%d")


class Outputs:

    # Rules for outputs
    # if output_root already exists, either 1) create a new ouput_root or 2) allow overwriting
    # all subdirectories will not be created if they already exists

    def __init__(self, output_root, overwrite=False):
        """ output_root : Full path of output directory. """

        if os.path.isdir(output_root):

            if not overwrite:
                suffix = 0
                while not os.path.isdir(output_root):
                    output_root = "{}_{}".format(output_root, suffix)
                    suffix += 1

                os.mkdir(output_root)
                print("directory already exists. new output directory created.")

        else:
            try:
                os.mkdir(output_root)
            except FileExistsError:
                pass

        self.root_dir = output_root
        self.paths = dict()

    def get(self, file_key):

        return self.paths[file_key]

    def add_output(
        self, output, output_path, mkdir=False, add_root=False, custom_root=None
    ):
        """update the dictionary of outputs
           * output -- name/label of the output
           * output_file -- full path of the output_file"""

        if add_root and not custom_root:
            full_path = os.path.join(self.root_dir, output_path)
        elif custom_root and not add_root:
            full_path = os.path.join(custom_root, output_path)
        elif custom_root and add_root:
            raise ValueError(
                "two root paths were provided when creating output files/dirs."
            )
        else:
            full_path = output_path

        self.paths[output] = full_path

        if mkdir:
            try:
                os.mkdir(full_path)
            except FileExistsError as e:
                pass


# def start_logger(log_file):
#
#     logger = logging.getLogger(__name__)
#     logger.setLevel(logging.DEBUG)
#
#     # format1 = logging.Formatter("%(levelname)s - %(asctime)s %(name)s line#:%(lineno)d --> %(message)s")
#     format1 = logging.Formatter("%(levelname)s - %(asctime)s - %(name)s line#:%(lineno)d --> %(message)s","%x %H:%M:%S")
#     format2 = logging.Formatter("%(message)s")
#
#     fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
#     fh.setLevel(logging.DEBUG)
#     fh.setFormatter(format1)
#
#     ch = logging.StreamHandler()
#     ch.setLevel(logging.INFO)
#     ch.setFormatter(format2)
#
#     logger.addHandler(fh)
#     logger.addHandler(ch)
#
#     logger.info(f"Running: {' '.join(sys.argv)}\n")
#     logger.info(f'Logging to {log_file}')
#
#     return logger


def safe_mkdir(path):

    if not os.path.isdir(path):
        os.mkdir(path)
    # else:
    #     print("Will overwrite since output already exists: {}".format(path))


def report_iferr(plinkcmd, decode_stdout, err_bool, logger=None, type="Error"):
    """ report error/warning if it exists  """

    if np.any(err_bool):

        if logger:
            logger.debug(
                f"{type} while running plink cmd:\n  {plinkcmd}\n  {decode_stdout[np.where(err_bool)[0][0]]}\n"
            )

        else:
            print(f"While running plink cmd:\n  {plinkcmd}")
            print(f"{decode_stdout[np.where(err_bool)[0][0]]}\n")


def decode_errors(plink_stdout, plinkcmd):

    decode_stdout = plink_stdout.decode("utf-8").splitlines()
    err_bool = [x.startswith("Error:") for x in decode_stdout]

    return decode_stdout, err_bool


def decode_warnings(plink_stdout, plinkcmd):

    decode_stdout = plink_stdout.decode("utf-8").splitlines()
    warn_bool = [x.startswith("Warning:") for x in decode_stdout]

    return decode_stdout, warn_bool


def error_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=None):

    # error handling
    decode_stdout, err_bool = decode_errors(plink_stdout, plinkcmd)
    decode_stderr, err_bool_stderr = decode_errors(plink_stderr, plinkcmd)

    # report
    report_iferr(plinkcmd, decode_stdout, err_bool, logger=logger, type="Error")
    report_iferr(plinkcmd, decode_stderr, err_bool_stderr, logger=logger, type="Error")

    any_errors = np.any(err_bool) | np.any(err_bool_stderr)

    return any_errors


def warning_check_plink_run(plink_stdout, plink_stderr, plinkcmd, logger=None):

    # get warnings
    warnings_stdout, warnings_bool = decode_warnings(plink_stdout, plinkcmd)
    warnings_stderr, warnings_bool_stderr = decode_warnings(plink_stderr, plinkcmd)

    # report
    report_iferr(
        plinkcmd, warnings_stdout, warnings_bool, logger=logger, type="Warning"
    )
    report_iferr(
        plinkcmd, warnings_stderr, warnings_bool_stderr, logger=logger, type="Warning"
    )

    any_warnings = np.any(warnings_bool) | np.any(warnings_bool_stderr)

    return any_warnings


def start_logger(log_file):

    logger = logging.getLogger("main")
    logger.setLevel(logging.DEBUG)

    # format1 = logging.Formatter("%(levelname)s - %(asctime)s %(name)s line#:%(lineno)d --> %(message)s")
    format1 = logging.Formatter(
        "%(levelname)s [%(asctime)s] - %(name)s line:%(lineno)d -> %(message)s",
        "%x %H:%M:%S",
    )
    format2 = logging.Formatter("%(message)s")

    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(format1)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(format2)

    logger.addHandler(fh)
    logger.addHandler(ch)

    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    logger.info(f"Logging to {log_file} start on {dt_string}")

    return logger


def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])


#
def report_mem():
    pid = os.getpid()
    process = psutil.Process(pid)

    rss = convert_size(np.round(process.memory_info().rss))
    # uss = convert_size(np.round(process.memory_info().uss))

    return f"PID:{pid} using {rss} rss"
