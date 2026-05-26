import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.cluster import HDBSCAN
from sklearn.cluster import DBSCAN
import argparse
import random

parser=argparse.ArgumentParser(description='Clusterize cmap using GREMLIN')
parser.add_argument('--matrix', required=True, help='Matrix of zscore')
parser.add_argument('--out', help='Output txt file to write the clusters (i.e., the residues to be masked in the MSA, one cluster per line)')
parser.add_argument('--plot', help='If you do not want to plot it and just save the output, write an input here')
parser.add_argument('--save_by_points', help='If you want to save the clusters by points, write the name of the file where you want to save it here')
args = parser.parse_args()

zscore = np.loadtxt(args.matrix, delimiter=' ')

def points(zscore, th, diag):
    data = []
    max_dist = len(zscore[0])/np.sqrt(2)
    for i in range(len(zscore[0])):
        for j in range(len(zscore[0])):
            dist = abs(i-j)/np.sqrt(2) # distancia a la diagonal
            th_weighted = (th-1)*(1-dist/max_dist)+1
            if zscore[i][j] > th_weighted and j not in range(i-diag, i+diag+1) and j > i:## if zscore[i][j] > 1.5 and j not in range(i-7, i+8) and j > i: ## zscore[i][j] > 2 and j not in range(i-3, i+4) and j > i
                data.append([i, j])
    return data


def clustering(data, cluster_size, epsilon):
    #data = np.array(data)
    #hdbscan = HDBSCAN(min_cluster_size=2, cluster_selection_epsilon=4)#(eps=4, min_samples=3), HDBSCAN(min_cluster_size=3, cluster_selection_epsilon=2)
    hdbscan = HDBSCAN(min_cluster_size=cluster_size, cluster_selection_epsilon=epsilon)
    labels = hdbscan.fit_predict(data)
    return labels


num_clusters = 0
for th in np.arange(1, 3, 0.05):
    data_new = points(zscore, th, 6)
    data_new = np.array(data_new)
    if len(data_new) < len(zscore[0])*0.10:
        break
    for size in range(2, 6):
        for eps in np.arange(2, 6, 0.2):
            labels_new = clustering(data_new, size, eps)
            unique_labels_new = set(labels_new)
            num_clusters_new = len(unique_labels_new) - (1 if -1 in labels_new else 0)  # Exclude noise (label == -1)
            #if num_clusters_new > num_clusters:
            if abs(num_clusters_new - 100) < abs(num_clusters - 100):
                change = True
                for lab in unique_labels_new:
                    if lab != -1:
                        i_values = [i for i, _ in data_new[labels_new == lab]]
                        j_values = [j for _, j in data_new[labels_new == lab]]
                        max_distance_i = max(i_values) - min(i_values)
                        max_distance_j = max(j_values) - min(j_values)

                        #if max_distance_i > len(zscore[0])*0.10 or max_distance_j > len(zscore[0])*0.10:
                        if max_distance_i > len(zscore[0])*0.15 or max_distance_j > len(zscore[0])*0.15:
                            change = False
                            break
                        
                if change == True:
                    labels = labels_new.copy()
                    num_clusters = num_clusters_new
                    unique_labels = unique_labels_new.copy()
                    data = data_new.copy()
                    best_eps = eps
                    best_size = size
                    best_th = th
            


#unique_labels = set(labels)
#num_clusters = len(unique_labels) - (1 if -1 in labels else 0)  # Exclude noise (label == -1)

print(f'Number of clusters: {num_clusters}\n')
print(f'Best th = {best_th}, best size = {best_size}, best eps = {best_eps}\n')

if args.out:
    original_clusters = []
    for lab in unique_labels:
        
        '''
        residues = []
        for i, point in enumerate(data):
            if labels[i] == lab:
                if point[0] not in residues:
                    residues.append(point[0])
                if point[1] not in residues:
                    residues.append(point[1])'
        '''
        ##### Ara farem que tots els residus que estiguin entre el màxim i el mínim d'un cluster siguin tots presos com a "coevolutius"
        if lab != -1:   # Excloem els punts de "soroll"
            residues_0 = []
            residues_1 = []
            for i, point in enumerate(data):
                if labels[i] == lab:
                    if point[0] not in residues_0:
                        residues_0.append(point[0])
                    if point[1] not in residues_1:
                        residues_1.append(point[1])
            max_0 = int(max(residues_0))
            min_0 = int(min(residues_0))
            max_1 = int(max(residues_1))
            min_1 = int(min(residues_1))
            residues_0 = list(range(min_0, max_0 + 1))
            residues_1 = list(range(min_1, max_1 + 1))
            residues = []
            for res_0 in residues_0:
                if res_0 not in residues:
                    residues.append(res_0)
            for res_1 in residues_1:
                if res_1 not in residues:
                    residues.append(res_1)
            ########## (és la part canviada per agafar tots els residus que hi ha entre el màxim i el mínim d'un cluster)

            #residues = [str(residue.item()) for residue in residues]
            residues = [str(residue) for residue in residues]
            original_clusters.append(residues)
            with open(args.out, 'a') as f:
                f.write(' '.join(residues) + '\n')

    '''
    # El següent for és per fer el merge dels clusters que es toquen entre ells
    # Merged clusters: we will use DBSCAN to see the clusters that appear so we can discard the merging that mask the majority of the sequence
    for k, clust in enumerate(original_clusters):
        new_cluster = clust.copy()
        total_dist = 0
        maximum = max(clust)
        minimum = min(clust)
        exclude_list = [k]
        too_long = False

        # Provem d'agafar els clusters que tenen overlap amb cada clsuter de forma random, i parem d'afegir-ne quan la longitud és massa gran (o quan ja hem fet tots els clusters)
        while len(exclude_list) < len(original_clusters) and too_long == False:
            valid_choices = [l for l in range(len(original_clusters)) if l not in exclude_list]
            random_number = random.choice(valid_choices)
            exclude_list.append(random_number)

            clust2 = original_clusters[random_number]
            for clust2 in original_clusters:
                if clust2 != clust:
                    interseccion = list(set(clust) & set(clust2))
                    if len(interseccion) > 0:
                        #new_cluster.extend(clust2)
                        union = set(new_cluster) | set(clust2)
                        new_cluster = list(union)


            db = DBSCAN(eps=10, min_samples=2)
            res = np.array(new_cluster, dtype = float)
            res = res.reshape(1, -1)
            #print(res)
            clusters = db.fit_predict(res)
            unique_clusters = set(clusters)
            for lab in unique_clusters:
                if lab != -1:
                    pts = res[clusters == lab]
                    minim = pts.min()
                    maxim = pts.max()
                    distancia = abs(maxim-minim)
                    total_dist = total_dist + distancia
                    if distancia > len(zscore[0])*0.17 or total_dist > len(zscore[0])*0.38:
                        too_long = True
                        break
        
        #print(new_cluster)
        new_cluster = [str(residue) for residue in new_cluster]
        with open(args.out, 'a') as f:
            f.write(' '.join(new_cluster) + '\n')
    '''
    



if not args.plot:
    colormap = plt.get_cmap('tab20', num_clusters)
    # List of marker styles to use for different clusters
    marker_styles = ['o', 's', '^', '*', 'D', 'v', '<', '>', 'P', 'H']

    # Plot the data with distinct colors and shapes for each cluster
    plt.figure(figsize=(6, 6))

    # Plot points belonging to clusters (label >= 0)
    for i in range(num_clusters):
        cluster_points = data[labels == i]
        plt.scatter(cluster_points[:, 0], cluster_points[:, 1], 
                    color=colormap(i), 
                    marker=marker_styles[i % len(marker_styles)],  # Cycle through marker shapes
                    label=f'Cluster {i}')

    # Plot noise points (label == -1) in a distinct color (black)
    noise_points = data[labels == -1]
    plt.scatter(noise_points[:, 0], noise_points[:, 1], 
                color='black', marker='x', label='Noise Points')

    # Customize the plot
    plt.title("HDBSCAN Clustering with Tab20 Colormap and Different Shapes")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.show()

# Output the cluster labels
print("Cluster labels:", labels)


if args.save_by_points:
    with open(args.save_by_points, 'w') as f:
        for label_clust_guardar in unique_labels:
            clust_guardar = []
            for i, point in enumerate(data):
                if labels[i] == label_clust_guardar:
                    clust_guardar.append(f'({point[0]},{point[1]})')
            f.write(','.join(clust_guardar) + '\n')
