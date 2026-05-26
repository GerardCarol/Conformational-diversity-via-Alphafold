import subprocess
import glob
import numpy as np
import re
import matplotlib.pyplot as plt
import argparse
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import NearestCentroid
from sklearn.metrics import pairwise_distances
import sys
from matplotlib.ticker import MaxNLocator
import os
import shutil
import json

parser=argparse.ArgumentParser(description='Program description')
parser.add_argument('--AF_dir', required=True, help='Directory of the AF generated structures')
parser.add_argument('--out_file', required=True, help='Output file to write the centroids')
parser.add_argument('--rmsd_dictionary', required=True, help='Matrix-like dictionary file (.json) of RMSD between proteins (the indices correspond to pdb file, i.e. 1REC.pdb)')
parser.add_argument('--out_dir', required=True, help='Output directory to write the cluster files')
args=parser.parse_args()

AF_dir = args.AF_dir
out = args.out_file
rmsd_file = args.rmsd_dictionary

with open(rmsd_file, "r") as archivo:
    rmsd_all_to_all = json.load(archivo)


matrix = [list(row.values()) for row in rmsd_all_to_all.values()] #Canvia de diccionari a matriu
matrix = np.array(matrix) #Si no es fa això després dóna errors perquè els índexs creus que no són enters
ind = list(rmsd_all_to_all.keys()) #Ens diu a quina posició (i.e., "i" i "j") correspon cada pdb (nomes ho fem un cop pq es una matriu quadrada,
# això és la bijecció entre index de fila o columna i pdb)

# cluster_th = 3.5
# canviar el th que vulguem, abans era 6
cluster_th = 6
clusters_agg = AgglomerativeClustering(n_clusters=None, distance_threshold=cluster_th, compute_full_tree=True, metric='precomputed', linkage='average').fit_predict(matrix) 

# Find the centroids
unique_labels = np.unique(clusters_agg)
'''
centroids = []
for label in unique_labels:
    cluster_indices = np.where(clusters_agg == label)[0]
    centroid = np.mean(matrix[cluster_indices], axis=0)
    centroids.append(centroid)

# Find the nearest protein to each centroid
nearest_proteins = []
for centroid in centroids:
    distances_to_centroid = np.linalg.norm(matrix - centroid, axis=1)
    nearest_idx = np.argmin(distances_to_centroid)
    nearest_protein = nearest_idx  # This should be your protein index in your dataset
    nearest_proteins.append(nearest_protein)'
'''

nearest_proteins = []
for label in unique_labels:
    # Get the indices of the proteins in the current cluster
    cluster_indices = np.where(clusters_agg == label)[0]
    
    # Calculate the average pairwise distance for each protein in the cluster
    cluster_distances = matrix[cluster_indices][:, cluster_indices]  # Submatrix of distances between proteins in the same cluster
    
    # The centroid is the average pairwise distance (mean of the submatrix)
    centroid = np.mean(cluster_distances, axis=1)  # Mean of each row (average distance to others in the cluster)

    # Step 3: Find the protein closest to the centroid of each cluster
    # To store the closest protein index for each cluster
    if len(centroid) == 1:
        nearest_proteins.append(cluster_indices[0])  # The only protein in the cluster is the nearest
        continue
    # Calculate the distance of each protein in the cluster to the centroid
    
    # Find the index of the protein closest to the centroid (smallest distance)
    nearest_protein_index = np.argmin(centroid)
    
    # Save the index of the nearest protein to the centroid
    nearest_proteins.append(cluster_indices[nearest_protein_index])


for protein in nearest_proteins:
    full_path = os.path.join(AF_dir, ind[protein])
    with open(out, 'a') as f:
        f.write(full_path + '\n')

with open(os.path.join(args.out_dir, 'clusters_agg.txt'), 'w') as c:
    for cluster in clusters_agg:
        c.write(f"{cluster}\n")


