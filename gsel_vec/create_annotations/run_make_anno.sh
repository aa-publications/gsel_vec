#!/bin/sh
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --array=1
#SBATCH --output=python_job_slurm_%A_%a.out
#SBATCH --job-name=dict




module load Anaconda3/4.4.0
source activate py36_r_ml
# anno_array=(dummy argweave  betascore  fst_eas_afr  fst_eur_afr  fst_eur_eas  gerp  iES_Sabeti  linsigh  phastCon100  phyloP100  xpehh_afr2_eas  xpehh_afr2_eur  xpehh_eas_eur)
anno_array=(dummy syn_age)

anno=${anno_array[${SLURM_ARRAY_TASK_ID}]}

echo $anno
python make_annotation_dictionary.py $anno