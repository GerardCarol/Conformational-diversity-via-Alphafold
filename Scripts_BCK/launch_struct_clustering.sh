#!/bin/bash

#SBATCH -e errors/err_struct/%j.err
#SBATCH -o outputs/out_struct/%j.out
#SBATCH --job-name struct
#SBATCH --time=10-00:00:00
#SBATCH -p mmb_cpu_sphr
#####SBATCH -p mmb_cpu_zen3
#SBATCH -c 1
#SBATCH --mem=2GB

module load Anaconda3
source activate AF2_TA_2

echo "Start at $(date)"

code=$1
####code="${TARGET_PDB}"
minus_code="${code,,}"

repo_dir=${PWD} 
dir_py=$repo_dir/Scripts/struct_msas

mkdir -p $repo_dir/struct_msas/proteins
dir=$repo_dir/struct_msas/proteins

for num in 1.5 2 2.5
do
	msa_dir=$repo_dir/gremlin/proteins/"$minus_code"_"$num"/clusters_new
	rmsd_file=$repo_dir/cluster/proteins/"$code"_"$num"/"$code"_rmsd.json

	cd $dir
	mkdir msas_"$minus_code"_"$num"

	python $dir_py/struct_clustering.py --cluster_th 6 --code $code --msa_dir $msa_dir --rmsd_file $rmsd_file --out_dir $dir/msas_"$minus_code"_"$num"
done
