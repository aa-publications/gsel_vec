# Evolutionary signatures from GWAS
This pipeline detects enrichment and depletion of evolutionary signatures from GWAS summary statistics.<br>
**Note: this pipeline is currently under development.**\

--> diagram place holder <---

## Install
**Please install these dependecies.** <br>1. [PLINK 1.90 beta](https://www.cog-genomics.org/plink2) <br>
2. [Anaconda](https://www.anaconda.com/products/individual)

**Run the** `commands` **as shown.** <br> 1) clone git repository: `git clone https://github.com/abraham-abin13/gsel_vec.git` <br>
2) change directory: `cd gsel_vec` <br>
4) run download_data.sh: `./download_data.sh` <br>
5) install package: `pip install .` <br>
6) create conda env: `conda env create -n conda-env -f conda_env.yml`

**Demo after installation** <br> 1) activate conda env `source activate manual_` <br> 
2) run pipeline with demo data. `top_gsel_vec.py`


## Run 
### Inputs 
The main command is `top_gsel_vec.py <analysis_name>   <gwas_summary_file>   <outputpath>`: 
- **analysis_name**: a string without spaces describing the analysis <br> 
- **gwas_summary_file**: full path to the gwas summary statistics file <br>
-  **outputpath**: full path to create output directory which will have the same name as **analysis_name**.

### Outputs
Outputs are saved to **outputpath** in a folder with the name **analysis_name**. The final outputs of the analyses are stored in the folder **final_outputs**. This contains... 




### Details
<gwas summary file format >
 **GWAS summary statistics file format:**<br>- should be a tab separated file with one row per SNP <br>- coordinates should be in GRChg37<br>- should have the following headers in this order: ['snp', 'chr', 'pos', 'p'] 
  



