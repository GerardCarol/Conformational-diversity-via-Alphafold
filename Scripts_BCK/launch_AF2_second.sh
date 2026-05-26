#!/bin/bash

#SBATCH -e errors/err_AF_second/%j.err
#SBATCH -o outputs/out_AF_second/%j.out
#SBATCH --job-name AF_2
#SBATCH --time=20-00:00:00
#SBATCH -p mmb_gpu_3090,irb_gpu_3090
#####SBATCH -p spot_gpu
#####SBATCH -p mmb_gpu_2080
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=20GB
################
#####SBATCH -p mmb_cpu_zen3
#####SBATCH -c 28

####code="${TARGET_PDB}"
#code=1DEG
code=$1
repo_dir=${PWD}

module load Anaconda3
source activate AF2_TA_2
export PYTHONPATH="/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/alphafold"
module load CUDA
	
mkdir -p $repo_dir/second_AF

for num in 1.5 2 2.5
do
	echo "Start at $(date)"

	minus_code="${code,,}"

	dir=/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/repositori
	msas_dir=$repo_dir/struct_msas/proteins/msas_"$minus_code"_"$num"
	outdir=$repo_dir/second_AF/AF_structures_"$minus_code"_"$num"

	num_elementos=$(ls -1 "$msas_dir" | wc -l)
	printf -v num_fin "%03d" $((num_elementos - 1))

	cd $repo_dir/second_AF
	mkdir -p AF_structures_"$minus_code"_"$num"

	files=""
	for i in $(seq -f "%03g" 0 $num_fin); do
	  files="$files $msas_dir/OUT_${code}_$i.fasta"
	done

	python $repo_dir/Scripts/AF2/RunAF2multi.py $files --model_num 1 --recycles 3 --output_dir $outdir --af2_dir /apps/alpha_fold

	echo "End at $(date)"
done
