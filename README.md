# Evolutionary signatures from GWAS
This pipeline detects enrichment and depletion of evolutionary signatures from GWAS summary statistics.<br>
**Note: this pipeline is currently under development.**

**Required Software:** 
1) [Python 3.6](https://www.python.org/downloads/release/python-360/)   
2) [PLINK 1.90 beta](https://www.cog-genomics.org/plink2)

**Required Data**   
1) Download both folders to the `gsel_vec` directory.  
a) "snpsnap_database"  - [link (2.4 GB)](https://drive.google.com/drive/folders/1P9r9axKakwY20eD_f3NCoRY0g1aLcp2T?usp=sharing)    
b) "1kg" - [link (1.5 GB)](https://drive.google.com/drive/folders/1yjp31LhZSi2Ftu_QmgKDKevLqHKJTH-0?usp=sharing)


2) Download "anno_dict" folder to `gsel_vec/create_annotations`.  
a) "anno_dict" [link (4.1 GB)](https://drive.google.com/drive/folders/1dps7iWshulKKEukxCdBu6MTy3j2s8KCj?usp=sharing)
  
3) Install conda environment to run scripts.  
	`conda env create -n conda-env -f evo_sig_gwas.yml`


# Inputs
**Pass these required arguments (in this order) to run_evo_sig.py**  <br>1) **analysis_name:** a string without spaces describing the name of the analysis <br> 2) **gwas_summary_file:** full path to gwas summary statistics file <br> 3) **outputpath:** full path to create output directory (with the same name as "analysis_name"

**GWAS summary statistics file format:**<br>- should be a tab separated file with one row per SNP <br>- coordinates should be in GRChg37<br>- should have the following headers in this order: ['snp', 'chr', 'pos', 'p']


## Example
See demo.sh for example.
