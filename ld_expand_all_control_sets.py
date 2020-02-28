#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import argparse
import time
import numpy as np
import pandas as pd
from datetime import datetime

from functools import partial
from glob import glob

import logging
logger = logging.getLogger('main.{}'.format(__name__))

sys.path.append("/scratch/abraha1/gsel_/gsel_vec")
from helper_general import Outputs

DATE = datetime.now().strftime('%Y-%m-%d')



# -----------
# FUNCTIONS
# -----------

def parse_input_args():


    # TODO: convert all args to requried args

    parser = argparse.ArgumentParser(description='Run expand control set.')

    parser.add_argument('lead_ld_counts_file',
                        action='store', type=str,
                        help="table of lead GWAS SNPs")

    parser.add_argument('gwas_snps_r2_file',
                        action='store', type=str,
                        help='table of r2 between gwas snps')

    parser.add_argument('matched_file',
                        action='store', type=str,
                        help='input/lead gwas snps with matched SNPs')

    parser.add_argument('control_ld_dir',
                        action='store', type=int,
                        help='directory with ld SNPs for control snps')

    parser.add_argument('output_dir',
                        action='store', type=str,
                        help="output_dir")

    # retrieve passed arguments
    args = parser.parse_args()
    lead_ld_counts_file = args.lead_ld_counts_file
    gwas_snps_r2_file = args.gwas_snps_r2_file
    matched_file = args.matched_file
    control_ld_dir = args.control_ld_dir
    output_dir = args.output_dir



    return lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, output_dir

def force_ld_control_snps(ld_df, controls_snps):
    # make sure that the control SNPs are all in 'SNP_A'


    # keep all control snps in "SNP_A"
    temp_a = ld_df.loc[ld_df['SNP_A'].isin(controls_snps)].copy()
    # keep all control snps in "SNP_B"
    temp_b = ld_df.loc[ld_df['SNP_B'].isin(controls_snps)].copy()
    # flip SNP A and SNP B assignments
    temp_b.columns = ['CHR_B', 'BP_B', 'SNP_B', 'CHR_A', 'BP_A', 'SNP_A', 'R2', 'DP']
    temp_b = temp_b.loc[:, temp_a.columns]


    forced_ld_df = pd.concat([temp_a, temp_b])

    logger.info("* Note: Out of {:,} control snps, {:,} ({:.2f}%) snps were not found after running plink ld.".format(len(controls_snps), len(controls_snps) - forced_ld_df.SNP_A.nunique(), (len(controls_snps) - forced_ld_df.SNP_A.nunique())/len(controls_snps) ))

    return forced_ld_df

def load_ld_df(LD_DIR):
    ld_df = pd.DataFrame()
    for ldfile in glob(LD_DIR+"/*.ld"):

        temp_ld_df = pd.read_csv( ldfile, sep="\s+")
        ld_df = ld_df.append(temp_ld_df)

    return ld_df

def shuffle_and_list_r2(x, max_n):
    "return a list with all ld snps for a given control snp. NOTE: the control snp is appended as the first element."
    # note: max_n is the length of the list to be return;
    #       this will be equal to the dimension of the mask array storing the LD SNPs, which should include the control SNP


    num_ld_snps = max_n-1 # since we will append the control snp

    control_snp =  [x.SNP_A.values[0]]
    ld_snps_of_control_snp =  x.SNP_B
    ldsnps_lst = list(ld_snps_of_control_snp)

    r2_of_ld_snps =  x.R2
    r2_lst =  list(r2_of_ld_snps)


    if len(ldsnps_lst) >= num_ld_snps:
        return [1] + r2_lst[0:num_ld_snps]      # I checked to make sure the incides & slices make sense

    else:

        mulpt = int(np.ceil(num_ld_snps/len(ldsnps_lst)))   # I checked to make sure that number of times it is tiled is correct

        return [1] + (r2_lst*mulpt)[0:num_ld_snps]

def shuffle_and_list(x, max_n):
    "return a list with all ld snps for a given control snp. NOTE: the control snp is appended as the first element."
    # note: max_n is the length of the list to be return;
    #       this will be equal to the dimension of the mask array storing the LD SNPs, which should include the control SNP


    num_ld_snps = max_n-1 # since we will append the control snp

    control_snp =  [x.SNP_A.values[0]]
    ld_snps_of_control_snp =  x.SNP_B
    ldsnps_lst = list(ld_snps_of_control_snp)

    # r2_of_ld_snps =  x.R2
    # r2_lst =  list(r2_of_ld_snps)


    if len(ldsnps_lst) >= num_ld_snps:
        return control_snp + ldsnps_lst[0:num_ld_snps]      # I checked to make sure the incides & slices make sense

    else:

        mulpt = int(np.ceil(num_ld_snps/len(ldsnps_lst)))   # I checked to make sure that number of times it is tiled is correct

        return control_snp + (ldsnps_lst*mulpt)[0:num_ld_snps]

def same_chrom_check(lead_control_ld_expanded_df):
    temp_chrom_df = pd.concat((lead_control_ld_expanded_df.loc[:, ['lead_snp','R2','lead_snp_bool']], lead_control_ld_expanded_df.loc[:, [f'Set_{x}' for x in np.arange(1,101)]].applymap(lambda x: x.split(":")[0])), axis=1)
    bool_df = temp_chrom_df.loc[:, temp_chrom_df.columns.difference(['R2','lead_snp_bool'])].replace('None',np.nan).groupby('lead_snp').nunique() == 1

    return bool_df, temp_chrom_df

def explode(df, lst_cols, fill_value='', preserve_index=False):

    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res

def to_row_per_csnp(selected_ld_snps, n_control_snps):
    raw_control_df = pd.DataFrame(selected_ld_snps.swapaxes(1,2).reshape(1,-1, n_control_snps).squeeze())
    raw_control_df.replace('', np.nan, inplace=True)
    control_df = raw_control_df[~raw_control_df.isnull().any(1)].copy()
    control_df.reset_index(inplace=True, drop=True)
    control_df.columns = ['Set_{}'.format(num) for num in np.arange(1, control_df.shape[1]+1)]
    return control_df

def ld_expand_all_control_snps(lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, output_root):

    tstart = time.time()

    # set up outputs
    logger.info("Starting to LD expand all control snps.")
    output_dir = os.path.join(output_root, 'ld_expanded_control_and_input_snps')
    OutObj = Outputs(output_dir, overwrite=True)
    OutObj.add_output('ld_expanded_output' ,'ld_expanded_all_control_sets.tsv', add_root=True)
    OutObj.add_output('ld_r2_expanded_output' ,'r2_ld_expanded_all_control_sets.tsv', add_root=True)

    # %%

    # lead_ld_counts_file='/scratch/abraha1/gsel_/gsel_vec/test/snp_list_output/testlist_clump/input_snps_with_binned_ldsnps_counts.tsv'
    # gwas_snps_r2_file='/scratch/abraha1/gsel_/gsel_vec/test/snp_list_output/testlist_clump/input_and_ld_snps.tsv'
    # matched_file='/scratch/abraha1/gsel_/gsel_vec/test/snp_list_output/testlist_matched_snps/matched_snps.tsv'
    # output_root="/scratch/abraha1/gsel_/gsel_vec/test/snp_list_output"
    # control_ld_dir='/scratch/abraha1/gsel_/gsel_vec/test/snp_list_output/testlist_get_ldsnps_for_control_snps/ld_snps_for_control_snps'


    # %%
    ###
    ###   set up ld thresholds
    ###

    # IMPORTANT: LD BINS MUST GO FROM HIGH TO LOW
    # ld_thresholds =['ld<=1.0', 'ld<=0.9']
    ld_thresholds =['ld<=1.0']



    ###
    ###   load matched control sets
    ###

    matched_df = pd.read_csv(matched_file, sep="\t")
    control_snps = matched_df.iloc[:, 1:].values.flatten()
    n_lead_snps = matched_df.shape[0]
    n_control_snps = matched_df.shape[1]-1

    logger.info("LD expanding control sets for {:,} lead snps with {:,} control snps.".format(n_lead_snps, n_control_snps))



    ###
    ###   load lead gwas snps and their LD snps; make a table
    ###

    # r2 of input gwas snps
    gwas_lead_ld_df = pd.read_csv(gwas_snps_r2_file, sep="\t")

    # number of LD snps required per lead gwas snps
    ld_counts_df = pd.read_csv(lead_ld_counts_file, sep="\t")
    lead_gwas_snps = ld_counts_df.lead_snp.unique().tolist() # determine the lead SNPs from the ld count table (since this includes input lead snps that don't have any snps in LD )


    # add row for lead SNP in the 'ld_snp' column (this is to get thr right format for the output)
    temp_lead_df = pd.DataFrame({'lead_snp': lead_gwas_snps, 'ld_snp': lead_gwas_snps, 'R2': [1.0]*len(lead_gwas_snps)})
    final_lead_ld_df = pd.concat((gwas_lead_ld_df, temp_lead_df),axis=0)
    final_lead_ld_df['lead_snp_bool'] = final_lead_ld_df['lead_snp'] == final_lead_ld_df['ld_snp']
    final_lead_ld_df.sort_values(['lead_snp','lead_snp_bool' ,'R2'] ,ascending=False, inplace=True)
    final_lead_ld_df.reset_index(inplace=True, drop=True)

    # define the order of lead snps
    ordered_lead_snps = final_lead_ld_df[~final_lead_ld_df.duplicated(subset=['lead_snp'], keep='first')].lead_snp.tolist()
    ordered_lead_snps_dict = dict(zip(ordered_lead_snps, np.arange(len(ordered_lead_snps))))

    # reorder the ld_counts_df rows (table of ld snps for input/gwas snp)
    temp_ld_counts_df = ld_counts_df.copy()
    temp_ld_counts_df.set_index('lead_snp', inplace=True)
    ordered_ld_counts_df = temp_ld_counts_df.reindex(ordered_lead_snps)




    ###
    ###   load ld pairs for the control snps and *shuffle*
    ###

    ld_df = load_ld_df(control_ld_dir)

    # force teh control snp to be in the first column
    forced_ld_df = force_ld_control_snps(ld_df, control_snps)
    forced_ld_df['ld_bin'] = 1/np.round(forced_ld_df.R2*10) # take the reciprocal so that it sorts in decreasing r2

    # SHUFFLE to simulate random selection later on...
    rand_ind = np.arange(forced_ld_df.shape[0])
    np.random.seed(22)
    np.random.shuffle(rand_ind)
    forced_ld_df['rand_int'] = rand_ind
    forced_ld_df.sort_values(['SNP_A','ld_bin','rand_int'], inplace=True)



    ###
    ###   create mask and ld snps 3-D tensor
    ###

    # create a mask that will select SNPs in LD with the control snps
    #       * inner most array will be a vector of 1 or 0, serving as the mask to select LD SNPs for a specific control snp
    #       * (the first element will always be 1 to include the control SNP)
    #  Example:
    #           # control snp 1        # control snp 2
    #   [   [   [1, 1, 0, 0, 1, 0],   [1, 1, 1, 1, 0, 0], ...   ],    --> input/lead gwas snp #1
    #       [   [1, 1, 0, 0, 1, 0],   [1, 1, 1, 1, 0, 0], ...   ], ...  ]  --> input/lead gwas snp #2



    ld_bins_expanded_df = pd.DataFrame()
    ld_bins_r2_expanded_df = pd.DataFrame()
    for ld_threshold in  ld_thresholds:
        # ld_threshold = ld_thresholds[0]


        uppper_ld_threshold = float(ld_threshold.split("=")[1])
        lower_ld_threshold = uppper_ld_threshold-0.1
        logger.info("\t matching ld threshold {} < r2 ≤ {}".format(lower_ld_threshold, uppper_ld_threshold))



        # -----------
        # create mask
        # -----------

        # create a binary vector for each input/lead snp where 1 is a placeholder to a required LD SNP
        inner_vec_length = ordered_ld_counts_df.loc[:, ld_threshold].max() + 1  # +1 so that the control snp is included as the first element
        padded_mask_per_input_snp = np.array( [  np.pad(np.ones(x+1), (0, inner_vec_length - (x+1)) , mode='constant') for x in ordered_ld_counts_df.loc[:, ld_threshold]])
        assert padded_mask_per_input_snp.shape[0] == n_lead_snps, "number of lead snps don't match with mask array "

        # tile equal to number of control sets and swap axes to ensure input snps get the outer axes
        ld_mask = np.swapaxes(np.tile(padded_mask_per_input_snp, (n_control_snps,1,1)), 0,1)
        assert ld_mask.shape[1]  == n_control_snps and ld_mask.shape[0]  == n_lead_snps, "array shapes don't match up with inputs"


        # -----------
        # create nested ld array
        # -----------

        # follow the same structure as teh nested_mask array
        # except inner most array will have LD SNPs instead of 1 or 0
        # note that each inner array should tile (this is to simulate sampling with replacement)
        # the length of the inner vector will be determined by the length of the innerrmost array of nested_mask


        # load ld pairs w/ *forced* the input/lead linked control_snp to be in teh first column
        # NOTE: * filter for r2 -->  we remove variants w/ r2 above the higher threhsold in the window
        #       * this means LD SNPs will be picked randomly from stratified LD bins that are decreasing in r2
        og_filt_ld_df = forced_ld_df.loc[ forced_ld_df['R2'] <= uppper_ld_threshold, ['SNP_A','SNP_B','R2']].copy()

        # remove aall SNP_A SNP_A r2=1 hits so that it doest not get tiled when there are not LD SNPs
        filt_r2_ld_df = og_filt_ld_df.loc[~og_filt_ld_df.duplicated(['SNP_A','SNP_B'], keep=False), :].copy()
        filt_ld_df = og_filt_ld_df.loc[~og_filt_ld_df.duplicated(['SNP_A','SNP_B'], keep=False), ['SNP_A','SNP_B']].copy()
        filt_ld_df.reset_index(drop=True, inplace=True)


        par_shuf = partial(shuffle_and_list, max_n = ld_mask.shape[2])
        par_r2_shuf = partial(shuffle_and_list_r2, max_n = ld_mask.shape[2])


        # for each input/lead snp
        all_snps = list()
        all_r2 = list()
        for input_snp in ordered_lead_snps:
            # input_snp = ordered_lead_snps[0]

            # get all control snps and remember the order!
            csnps_for_isnp = matched_df.loc[matched_df['lead_snp'] == input_snp, ['Set_{}'.format(x) for x in range(1,n_control_snps+1) ]].values.flatten()

            uniq_csnps = np.unique(csnps_for_isnp)
            csnps_indices =  dict(zip(uniq_csnps, [np.where(csnps_for_isnp == uniq_csnp)[0].tolist() for uniq_csnp in uniq_csnps]))


            # update r2 table with control snps that are not in it with 'None' as a placeholder
            missing_control_snps =  set(csnps_for_isnp) - set(filt_ld_df.SNP_A.unique())
            miss_df = pd.DataFrame({'SNP_A': list(missing_control_snps), 'SNP_B':'None'})
            miss_r2_df = pd.DataFrame({'SNP_A': list(missing_control_snps), 'SNP_B':'None', 'R2':'None'})

            updated_ld_df = pd.concat((miss_df, filt_ld_df))
            updated_r2_ld_df = pd.concat((miss_r2_df, filt_r2_ld_df))

            # collect all the r2 for LD snps into one row
            ld_r2_array_df = updated_r2_ld_df.loc[ updated_r2_ld_df['SNP_A'].isin(csnps_for_isnp)].groupby('SNP_A', sort=False).apply(par_r2_shuf).reset_index()
            ld_r2_array_df.rename({0:'r2_ld_array'}, axis=1, inplace=True)

            # collect all the D snps into one row
            ld_array_df = updated_ld_df.loc[ updated_ld_df['SNP_A'].isin(csnps_for_isnp)].groupby('SNP_A', sort=False).apply(par_shuf).reset_index()
            ld_array_df.rename({0:'ld_array'}, axis=1, inplace=True)

            # order control snps, this is necessary to preserve order between different ld bin iterations
            ld_array_df['order'] = ld_array_df.SNP_A.map(csnps_indices)
            ordered_ld_array_df = explode(ld_array_df, 'order')
            ordered_ld_array_df.sort_values('order', inplace=True)

            # order for the r2 dataframe as well
            ld_r2_array_df['order'] = ld_r2_array_df.SNP_A.map(csnps_indices)
            ordered_ld_r2_array_df = explode(ld_r2_array_df, 'order')
            ordered_ld_r2_array_df.sort_values('order', inplace=True)


            # extract arrays of ld snps for each control snp
            all_snps.append(ordered_ld_array_df.loc[:,'ld_array'].values.tolist())
            all_r2.append(ordered_ld_r2_array_df.loc[:,'r2_ld_array'].values.tolist())



        # select LD SNPs based on ld_mask
        ld_snps_array = np.array(all_snps, dtype='object', order='F')
        ld_snps_r2_array = np.array(all_r2, dtype='object', order='F')
        selected_ld_snps = np.where( np.array(ld_mask, dtype=bool), ld_snps_array, '')
        selected_ld_r2_snps = np.where( np.array(ld_mask, dtype=bool), ld_snps_r2_array, '')


        # from one column per control set, to ...
        # one row per control snp (lead control snp being first row in each group)
        control_df = to_row_per_csnp(selected_ld_snps,n_control_snps)
        control_r2_df = to_row_per_csnp(selected_ld_r2_snps, n_control_snps)


        # create lead snp columns to concat w/ control snps df
        r2_filt_final_lead_ld_df = final_lead_ld_df.loc[  (final_lead_ld_df['R2'] > lower_ld_threshold)  & (final_lead_ld_df['R2'] <= uppper_ld_threshold)].reset_index(drop=True)


        assert control_df.shape[0] == r2_filt_final_lead_ld_df.shape[0], "Number of rows do not match between the two dataframes (lead snps and control snps) that are going to be concatentated "
        assert control_r2_df.shape[0] == r2_filt_final_lead_ld_df.shape[0], "Number of rows do not match between the two dataframes (lead snps and control snps) that are going to be concatentated "


        if (uppper_ld_threshold != 1):
            # add back the row containg the lead SNP
            # need to do this since,  3-D tensor's first element is the control snps for which the LD SNPs belong to
            temp_lead_snps_df = final_lead_ld_df.loc[ final_lead_ld_df['lead_snp_bool']==True, :].copy()
            r2_filt_final_lead_ld_df= r2_filt_final_lead_ld_df.append(temp_lead_snps_df)

            # ensure that the order remains teh same; some SNPs may have fell through if they didn't meet the r2 LD SNPS criteria.
            r2_filt_final_lead_ld_df['temp_order'] = r2_filt_final_lead_ld_df.lead_snp.map(ordered_lead_snps_dict)*-1
            r2_filt_final_lead_ld_df.sort_values(['temp_order','lead_snp_bool','R2'], ascending=False, inplace=True)
            r2_filt_final_lead_ld_df.reset_index(inplace=True, drop=True)
            r2_filt_final_lead_ld_df.drop('temp_order', axis=1, inplace=True)



        # concat lead snp and control snps data frames
        lead_control_ld_expanded_df = pd.concat( (r2_filt_final_lead_ld_df, control_df), axis=1, ignore_index=False)
        lead_control_ld_expanded_r2_df = pd.concat( (r2_filt_final_lead_ld_df, control_r2_df), axis=1, ignore_index=False)
        ld_bins_expanded_df = ld_bins_expanded_df.append(lead_control_ld_expanded_df)
        ld_bins_r2_expanded_df = ld_bins_r2_expanded_df.append(lead_control_ld_expanded_r2_df)


    # remove the extra rows used to keep the control snps in order...
    dedup_df = ld_bins_expanded_df[~ld_bins_expanded_df.duplicated(subset=['lead_snp','ld_snp','R2'], keep='first')].sort_values(['lead_snp','ld_snp'])
    dedup_df.sort_values(['lead_snp','lead_snp_bool'], ascending=False, inplace=True)

    dedup_r2_df = ld_bins_r2_expanded_df[~ld_bins_r2_expanded_df.duplicated(subset=['lead_snp','ld_snp','R2'], keep='first')].sort_values(['lead_snp','ld_snp'])
    dedup_r2_df.sort_values(['lead_snp','lead_snp_bool'], ascending=False, inplace=True)

    # round
    dedup_r2_df.replace('None', np.nan, inplace=True)
    dedup_r2_df = dedup_r2_df.round(2)
    dedup_r2_df.fillna('None', inplace=True)


    # write
    dedup_df.to_csv(OutObj.get('ld_expanded_output'), sep="\t", index=False)
    dedup_r2_df.to_csv(OutObj.get('ld_r2_expanded_output'), sep="\t", index=False)
    logger.info("Wrote expanded control snps to {}".format(OutObj.get('ld_expanded_output')))
    logger.info("Done LD expanding control snps in {:.2f} minutes".format((time.time()-tstart)/60))

    return OutObj

# -----------
# MAIN
# -----------


if __name__ == '__main__':


    lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, output_root = parse_input_args()
    ld_expand_all_control_snps(lead_ld_counts_file, gwas_snps_r2_file, matched_file, control_ld_dir, output_root)
