# Evolutionary signatures from GWAS
This pipeline will detect evolutionary signatures from GWAS summary statistics.

**Dependencies:**

1) Download the following folders to gsel_vec directory.
    - "snpsnap_database" - [link (2.4 GB)](https://drive.google.com/drive/folders/1P9r9axKakwY20eD_f3NCoRY0g1aLcp2T?usp=sharing)
    - "1kg" - [link (1.5 GB)](https://drive.google.com/drive/folders/1yjp31LhZSi2Ftu_QmgKDKevLqHKJTH-0?usp=sharing)

2) Download and move "anno_dict" folder to gsel_vec/create_annotations. 
    - "anno_dict" [link (4.1 GB)](https://drive.google.com/drive/folders/1dps7iWshulKKEukxCdBu6MTy3j2s8KCj?usp=sharing)
   

3) Install conda environment to run scripts.

`conda env create -n conda-env -f evo_sig_gwas.yml`


## Inputs
**Pass these required arguments (in this order) to run_evo_sig.py**
1) **analysis_name:** a string without spaces describing the name of the analysis
2) **gwas_summary_file:** full path to gwas summary statistics file
3) **outputpath:** Outputs will be stored in a directory (with the same nams as "analysis_name") under this path.



**GWAS summary statistics file format:**
- should be a tab separated file with one row per SNP
- coordinates should be in GRChg37
- should have the following headers in this order: ['snp', 'chr', 'pos', 'p']


## Example
See demo.sh for example.
