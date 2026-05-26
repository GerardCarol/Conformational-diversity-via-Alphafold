import argparse
import os
import glob
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import json

parser=argparse.ArgumentParser(description='Program description')

parser.add_argument('--cluster_th', required=True, help='cluster threshold (rmsd)')
parser.add_argument('--code', required=True, help='Original pdb that corresponds, i.g. 1DEG, 1JSA, etc.')
parser.add_argument('--out_dir', required=True, help='Output directory to write the cluster files')
parser.add_argument('--msa_dir', required=True, help='Directory where the msa files are located')
parser.add_argument('--rmsd_file', required=True, help='RMSD file to be used for clustering')
args=parser.parse_args()

code = args.code
minus_code = code.lower()
cluster_th=float(args.cluster_th)
out_dir = args.out_dir
rmsd_file = args.rmsd_file
msa_dir = args.msa_dir


with open(rmsd_file, "r") as archivo:
    rmsd_all_to_all = json.load(archivo)


matrix = [list(row.values()) for row in rmsd_all_to_all.values()] #Canvia de diccionari a matriu
matrix = np.array(matrix) #Si no es fa això després dóna errors perquè els índexs creus que no són enters
ind = list(rmsd_all_to_all.keys()) #Ens diu a quina posició (i.e., "i" i "j") correspon cada pdb (nomes ho fem un cop pq es una matriu quadrada,
# això és la bijecció entre index de fila o columna i pdb)

clusters_agg = AgglomerativeClustering(n_clusters=None, distance_threshold=cluster_th, compute_full_tree=True, metric='precomputed', linkage='average').fit_predict(matrix)
print(clusters_agg)

'''
with open(os.path.join(out_dir, 'clusters_agg.txt'), 'w') as c:
    for cluster in clusters_agg:
        c.write(f"{cluster}\n")
'''
        
unique_labels = np.unique(clusters_agg)
'''
for label in unique_labels:
    msas = []
    for i in range(0, len(clusters_agg)):
        if clusters_agg[i] == label:
            msas.append(i)
    
    ind = f"{msas[0]:03}"
    indj = f"{label:03}"
    with open(f'{msa_dir}/OUT_{code}_{ind}.a3m', 'r') as msa:
        lines = msa.readlines()
    with open(f'{out_dir}/OUT_{code}_{indj}.a3m', 'a') as out:    
        out.writelines(lines)

    for i in range(1, len(msas)):
        ind = f"{msas[i]:03}"
        with open(f'{msa_dir}/OUT_{code}_{ind}.a3m', 'r') as msa:
            lines = msa.readlines()
        sequences = lines[2:]
        with open(f'{out_dir}/OUT_{code}_{indj}.a3m', 'a') as out:
            out.writelines(sequences)
'''
for label in unique_labels:
    msas = []
    for i in range(len(clusters_agg)):
        if clusters_agg[i] == label:
            msas.append(i)
    
    indj = f"{label:03}"
    out_file = os.path.join(out_dir, f'OUT_{code}_{indj}.fasta')
    first, _ = os.path.splitext(ind[msas[0]])
    read_file = os.path.join(msa_dir, f'{first}.fasta')
    with open(read_file, 'r') as read:
        lines = read.readlines()
    with open(out_file, 'a') as out:    
        out.writelines(lines)
    
    for i in range(1, len(msas)):
        fastas, _ = os.path.splitext(ind[msas[i]])
        read_fastas = os.path.join(msa_dir, f'{fastas}.fasta')
        with open(read_fastas, 'r') as read:
            lines = read.readlines()
        sequences = lines[2:]
        with open(out_file, 'a') as out:
            out.writelines(sequences)