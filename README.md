# Evolutionary signatures from GWAS
This pipeline detects enrichment and depletion of evolutionary signatures from GWAS summary statistics.<br>
**Note: this pipeline is currently under development.**


![](<https://github.com/abraham-abin13/gsel_vec/blob/master/schematic/for_repo.jpg>)

## Install
**Please install these dependecies.**
<br>1. [PLINK 1.90 beta](https://www.cog-genomics.org/plink2)<br>2. [Anaconda](https://www.anaconda.com/products/individual) <br>
- _note: the shell should have path access to run the plink command_ 

**Run the** `commands` **as shown.** <br>
1) clone git repository: `git clone https://github.com/abraham-abin13/gsel_vec.git` <br>
2) change directory: `cd gsel_vec` <br>
3) *run download_data.sh: `./download_data.sh` <br>
- _note: requires the the 'wget' utility; if on a mac, install wget using brew (`brew install wget`)_

4) create conda env: `conda env create -f conda_env_manual_dask.yml` <br>
5) activate conda env: `source activate evo_dask` <br>
5) install package: `pip install .` <br>


## Test Suite After Installation <br>
1) activate conda env `source activate evo_dask` <br>
2) run the command `top_gsel_vec.py`
3) outputs are saved to a folder called 'testing' in the current directory. 
4) See Outputs section and Details sections below for interpreting output files. 


## Run
### Inputs
The main command is `top_gsel_vec.py <analysis_name>   <gwas_summary_file>   <output_path>`:
- **analysis_name**: a label describing the analysis <br>
- **gwas_summary_file**: full path to the gwas summary statistics file <br>
-  **output_path**: full path to create output directory which will have the same name as **analysis_name**.

### Outputs
Outputs are saved to **output_path** in a folder with the name **analysis_name**. The final outputs of the analyses are stored in the folder **final_outputs** with the following folders:
- **trait_enrichments/trait_enrichment.tsv**:  enrichment for each evolutionary measure at the trait level
-**pval_zscore_per_anno/*.tsv**: Z-score and empirical p-value for each trait associated locus compared to matched background; one file per evolutionary measure
- **plot_outputs/**: [to do]



### Details
A. GWAS Summary Statistics File Format
- should be a tab separated file with one row per SNP
- coordinates should be in GRChg37
- should have the following headers in this order: ['snp', 'chr', 'pos', 'p']

B. Interpreting **trait_enrichment.tsv** output
- one row per evolutionary measure

| Column  | Description |
| ------------- | ------------- |
| annotation  | name of evolutionary measure  |
| mean_trait_loci  | mean across all trait associated loci  |
| median_trait_loci  | median across all trait associated loci  |
| emp_pvalue  |  p-value from comparing to matched background   |
| \*_matched_*  |  matched distributions mean, median, std, 5th, & 95th percentiles   |
| n_lead_loci  |  number of trait associated regions after LD clumping   |
| enrich_per_mean_diff_by_genomstd  |  enrichment calculated as described below   |

* For each evolutionary measure, enrichment is calculated by subtracting the mean across trait-associated loci from mean across matched loci and dividing this by the genome-wide standard deviation for that evolutoinary measure:<br> <br><img src="https://render.githubusercontent.com/render/math?math=enrichment = \frac{trait_{mean} - matched_{mean}}{Genome\:wide_{std}}">

C. Interpreting **pval_zscore_per_anno/*.tsv** output
- one file per evolutionary measure
- one row per LD independent trait-assocaited region

| Column  | Description |
| ------------- | ------------- |
| lead_snp  | SNP with lowest P-value for a trait-associated region  |
| z_score  | z-score of trait-associated region compared to matched background  |
| lead_snp_anno  | mean of evolutionary measure of SNPs in a trait-associated region  |
| mean_controls  | mean of evolutionary measure of matched background  |
| pvalue  | empirical p-value from comparing evolutinoary measure at a trait-associated region (lead_snp_anno) to matched background   |
| corrected_pval_benj_hoch  | empirical p-value corrected for number of trait-associated regions  |
| annotation  | name of evoluationary measure   |
| final_qc_include  | if False, remove trait-associated region from downstream analysis  |

D. Methods
- For details please see manuscript.

### Running one chromosome at a time
For GWAS with a large number of trait-associated regions, this pipeline can be run on individual chromosomes. The results from each chromosome can be combined using the script `gsel_vec/gsel_vec/scripts/by_chr_calc_enrich.py`.
* to run this script, you need to have the dask python package installed; you can install the `conda_env_dask.yml` conda enviornment for this.


A. Instructions
1. split summary statistics into files for each chromosome
2. per chromosome summary stats files should be named chr<#>_<analysis_name>
3. create a folder for output_path named <analysis_name>
4. run `top_gsel_vec.py -a <analysis_name> -g <path_to_per_chromosome_summary_stats_file> -o <output_path> --run_by_chr` for each chromosome
    * note the `--run_by_chr` flag is added
    * note the <output_path> should be the full path to the folder named <analysis_name> created in step 3.
    * repeat this command with each chromosome file (changing the -g flag input) 
5. run `top_gsel_vec.py -a <analysis_name> -o <output_path> --run_by_chr_enrich`.
    * This will calculate trait wide enrichments.  
7. Trait-wide enrichment will saved to `<analysis_name>_extreme_regions_mean_enrichment_all_annotation.tsv`
