# Evolutionary signatures from GWAS
This pipeline will detect evolutionary signatures from GWAS summary statistics.



## Inputs
Pass these required arguments (in this order) to run_evo_sig.py
1) analysis_name: name of the analysis
2) gwas_summary_file: full path to gwas summary statistics file
3) outputpath: Outputs will be stored in a directory (with the same nams as "analysis_name") under this path.

GWAS Summary Stat Format:
- should be a tab separated file with one row per SNP
- should have the following headers in this order: ['snp','chr','pos','p']
