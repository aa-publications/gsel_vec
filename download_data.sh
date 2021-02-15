#!/bin/sh


# dowload required data and put them into the appropriate directories

ANNO_DICT="https://ndownloader.figshare.com/articles/13664564/versions/1"
SNPSNAP="https://ndownloader.figshare.com/files/26226920"
KG_DICT="https://ndownloader.figshare.com/articles/13664570/versions/1"

# OUTPUT_DIR="/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/gsel_vec/data"
OUTPUT_DIR="/dors/capra_lab/users/abraha1/temp_/test_dl_data"



wget -O ${OUTPUT_DIR}"/anno_dict/anno_dict.zip" $ANNO_DICT
unzip ${OUTPUT_DIR}"/anno_dict/anno_dict.zip" -d ${OUTPUT_DIR}/"anno_dict"
rm ${OUTPUT_DIR}"/anno_dict/anno_dict.zip"



output_name=${OUTPUT_DIR}"/snpsnap_database/snpsnap_database_ld0.9.tab.gz"
wget -O $output_name $SNPSNAP -d ${OUTPUT_DIR}/"snpsnap_database"


output_name=${OUTPUT_DIR}"/1kg/1kg.zip"
wget -O $output_name $KG_DICT
unzip $output_name -d ${OUTPUT_DIR}/"1kg"
rm $output_name


