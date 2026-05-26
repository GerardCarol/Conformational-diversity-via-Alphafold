#!/bin/bash

#SBATCH -e errors/err_cluster/%j.err
#SBATCH -o outputs/out_cluster/%j.out
#SBATCH --job-name clustering
#SBATCH --time=3-00:00:00
#####SBATCH -p mmb_cpu_zen3
#SBATCH -p mmb_cpu_sphr
#SBATCH -c 3
#SBATCH --mem=8GB

module load Anaconda3
source activate AF2_TA_2

echo "Start at $(date)"

PDB=$1
####PDB="${TARGET_PDB}"
pdb="${PDB,,}"

repo_dir=${PWD}

mkdir -p $repo_dir/cluster/proteins

for num in 1.5 2 2.5
do
	dir_af=$repo_dir/first_AF/AF_structures_"$pdb"_"$num"
	#dir=$repo_dir/cluster/scripts_cluster
	out_dir=$repo_dir/cluster/proteins

	cd $out_dir
	mkdir -p "$PDB"_"$num"

	python $repo_dir/Scripts/cluster/rmsd_for_cluster.py --dir $dir_af --rmsd_file $out_dir/"$PDB"_"$num"/"$PDB"_rmsd.json

	python $repo_dir/Scripts/cluster/cluster.py --AF_dir $dir_af --out_file $out_dir/"$PDB"_"$num"/"$pdb"_analysis_centroids_list.txt --rmsd_dictionary $out_dir/"$PDB"_"$num"/"$PDB"_rmsd.json --out_dir $out_dir/"$PDB"_"$num"
done
