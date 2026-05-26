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

unique_labels = np.unique(clusters_agg)

for label in unique_labels:
    msas = []
    for i in range(len(clusters_agg)):
        if clusters_agg[i] == label:
            msas.append(i)
    
    indj = f"{label:03}"
    out_file = os.path.join(out_dir, f'OUT_{code}_{indj}.fasta')
    first, _ = os.path.splitext(ind[msas[0]])
    read_file = os.path.join(msa_dir, f'{first}.fasta')

    total_msa = []
    for msa in msas:
        fastas, _ = os.path.splitext(ind[msa])

        ### Canvi fet per si no va bé  en majúscula o minúscula
        parts_fastas = fastas.split('_')
        read_fastas = os.path.join(msa_dir, f'{parts_fastas[0]}_{parts_fastas[1].upper()}_{parts_fastas[2]}.fasta')
        #####read_fastas = os.path.join(msa_dir, f'{fastas}.fasta')
        if not os.path.exists(read_fastas):
            read_fastas = os.path.join(msa_dir, f'{parts_fastas[0]}_{parts_fastas[1].lower()}_{parts_fastas[2]}.fasta')

        with open(read_fastas, 'r') as read:
            lines = read.readlines()
            lines = [line.strip() for line in lines]
        for i, line in enumerate(lines):
            if line.startswith('>') and line not in total_msa:
                total_msa.append(line)
                total_msa.append(lines[i+1])
    
    write_msa = [line + '\n' for line in total_msa]
    with open(out_file, 'w') as out:
        out.writelines(write_msa)
