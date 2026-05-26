#!/bin/bash

#SBATCH -e err_get_msas/%j.err
#SBATCH -o out_get_msas/%j.out
#SBATCH --job-name msas
#SBATCH --time=01:00:00
#SBATCH -p mmb_cpu_zen3
#####SBATCH -p mmb_cpu_sphr
#SBATCH --mem=2GB
#SBATCH -c 1

#source activate alphaflow
module load Anaconda3
#source activate AlphaFlow

echo "Start at $(date)"

### Per canviar de proteïna, només canviar el nom, no els directoris
pdb=$1

#/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/alphaflow/scripts.mmseqs_query
python Scripts/mmseqs_query.py --split $pdb.csv --outdir msas

echo "End at $(date)"

