#!/bin/bash

#SBATCH -e errors/err_gremlin/%j.err
#SBATCH -o outputs/out_gremlin/%j.out
#SBATCH --job-name coev
#SBATCH --time=10-00:00:00
#####SBATCH -p mmb_cpu_zen3
#SBATCH -p mmb_cpu_sphr
#SBATCH -c 8
#SBATCH --mem=32GB

module load Anaconda3
source activate AF2_TA

echo "Start at $(date)"

####ID="${TARGET_PDB}"
ID=$1
id="${ID,,}"

echo "$ID"

#dir=/data/mmb/TransAtlas/Gerard/MontseScripts/transatlas_af/prova/repositori
dir=${PWD}

cd $dir/msas
cd $(find . -type d -name "${id}*" -print -quit)
cd a3m
msa=$(realpath $(find . -type f -name "${id}*" -print -quit))

mkdir -p $dir/gremlin/proteins
for th in 1.5 2 2.5
do
	cd $dir/gremlin/proteins
	mkdir "$id"_"$th"
done

cd $dir/gremlin/proteins/""$id"_1.5"

linea=$(sed -n '2p' $msa)
longitud=${#linea}
python $dir/Scripts/gremlin/len.py --msa $msa --out_si si.a3m --out_no no.a3m --len $longitud

#source activate muscle
mafft --addfull no.a3m --keeplength --auto si.a3m > prova.fasta
#conda deactivate

#source activate gremlin
$dir/Scripts/gremlin/gremlin_cpp -i prova.fasta -o prova.txt -gap_cutoff 1 -mrf_o mrf.txt
#conda deactivate

#source activate plot
python $dir/Scripts/gremlin/plot_coev.py --matrix prova.txt --out zscore.txt --plot no
#conda deactivate


cd $dir/gremlin/proteins/""$id"_1.5"
#source activate analysis_filter
python $dir/Scripts/gremlin/cluster_cmap.py --matrix zscore.txt --out new_clusters.txt --plot no
mkdir clusters_new
#conda deactivate

#source activate plot
python $dir/Scripts/gremlin/cluster_mrf.py --mrf mrf.txt --msa prova.fasta --clusters new_clusters.txt --code $ID --out_dir clusters_new --th 1.5
#conda deactivate

for th in 2 2.5
do
	grem_dir=$dir/gremlin/proteins/""$id"_1.5"
	cp $grem_dir/prova.fasta $dir/gremlin/proteins/"$id"_"$th"
	cp $grem_dir/zscore.txt $dir/gremlin/proteins/"$id"_"$th"
	cp $grem_dir/mrf.txt $dir/gremlin/proteins/"$id"_"$th"
	cp $grem_dir/new_clusters.txt $dir/gremlin/proteins/"$id"_"$th"

	cd $dir/gremlin/proteins/"$id"_"$th"
	mkdir clusters_new

	#source activate plot
	python $dir/Scripts/gremlin/cluster_mrf.py --mrf mrf.txt --msa prova.fasta --clusters new_clusters.txt --code $ID --out_dir clusters_new --th $th
	#conda deactivate
done

echo "End at $(date)"
