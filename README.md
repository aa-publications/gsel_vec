# Evolutionary signatures from GWAS
This pipeline will detect evolutionary signatures from GWAS summary statistics.

**Dependencies:**

Please download the following folders to gsel_vec directory.
1) "snpsnap_database" - link (X GB)
2) "1kg" - link (X GB)
3) "create_annotations" - link (X GB)

Please use the conda enviornment provided.

`conda env create -n conda-env -f evo_sig_gwas.yml`


## Inputs
Pass these required arguments (in this order) to run_evo_sig.py
1) analysis_name: a string without spaces describing the name of the analysis
2) gwas_summary_file: full path to gwas summary statistics file
3) outputpath: Outputs will be stored in a directory (with the same nams as "analysis_name") under this path.



** GWAS summary statistics file format: **
- should be a tab separated file with one row per SNP
- should have the following headers in this order: ['snp','chr','pos','p']

