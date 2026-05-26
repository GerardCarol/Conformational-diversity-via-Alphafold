#!/bin/bash

#SBATCH -e errors/err_AF_first/%j.err
#SBATCH -o outputs/out_AF_first/%j.out
#SBATCH --job-name AF
#SBATCH --time=7-00:00:00
#SBATCH -p mmb_gpu_3090,irb_gpu_3090
#####SBATCH -p spot_gpu
#####SBATCH -p irb_gpu_3090
#####SBATCH -p mmb_gpu_2080
#SBATCH --gres=gpu:1
#SBATCH -c 3
#SBATCH --mem=20GB
################
#####SBATCH -p mmb_cpu_zen3
#####SBATCH -c 28


####code="${TARGET_PDB}"
#code=1DEG
code=$1
repo_dir=${PWD}

module load Anaconda3
#source activate AF2_TA
#source activate AF2
source activate AF2_TA_2
export PYTHONPATH="/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/alphafold"
module load CUDA

for num in 1.5 2 2.5
do

	echo "Start at $(date)"

	minus_code="${code,,}"
	
	#repo_dir="/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/repositori"
	dir=$repo_dir/gremlin/proteins/"$minus_code"_"$num"
	msas_dir=$dir/clusters_new
	mkdir -p $repo_dir/first_AF
	outdir=$repo_dir/first_AF

	cantidad_lineas=$(wc -l < "$dir/new_clusters.txt")
	printf -v num_fin "%03d" $((cantidad_lineas - 1))

	cd $outdir
	mkdir -p AF_structures_"$minus_code"_"$num"

	files=""
	for i in $(seq -f "%03g" 0 $num_fin); do
		files="$files $msas_dir/OUT_${code}_$i.fasta"
	done

	python $repo_dir/Scripts/AF2/RunAF2multi.py $files --model_num 1 --recycles 3 --output_dir $outdir/AF_structures_"$minus_code"_"$num" --af2_dir /apps/alpha_fold

	echo "End at $(date)"
done
