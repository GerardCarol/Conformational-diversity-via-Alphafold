import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import HDBSCAN
import argparse
import random


parser=argparse.ArgumentParser(description='Clusterize cmap using GREMLIN')
parser.add_argument('--matrix', required=True, help='Matrix of zscore')
parser.add_argument('--out', help='Output txt file to write the clusters (i.e., the residues to be masked in the MSA, one cluster per line)')
parser.add_argument('--plot', help='If you do not want to plot it and just save the output, write an input here')
args = parser.parse_args()

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
    hdbscan = HDBSCAN(min_cluster_size=cluster_size, cluster_selection_epsilon=epsilon)
    labels = hdbscan.fit_predict(data)
    return labels


if __name__ == "__main__":
    zscore = np.loadtxt(args.matrix, delimiter=' ')

    num_clusters = 0
    for th in np.arange(1, 3, 0.05):
        data_new = points(zscore, th, 6)
        data_new = np.array(data_new)
        if len(data_new) < len(zscore[0])*0.10: # This value can be changed to avoid too few points
            break
        for size in range(2, 6):
            for eps in np.arange(2, 6, 0.2):
                labels_new = clustering(data_new, size, eps)
                unique_labels_new = set(labels_new)
                num_clusters_new = len(unique_labels_new) - (1 if -1 in labels_new else 0)  # Exclude noise (label == -1)
                #if num_clusters_new > num_clusters:
                if abs(num_clusters_new - 100) < abs(num_clusters - 100):
                    labels = labels_new.copy()
                    num_clusters = num_clusters_new
                    unique_labels = unique_labels_new.copy()
                    data = data_new.copy()
                    best_eps = eps
                    best_size = size
                    best_th = th
    

    if len(unique_labels) == 1 and -1 in unique_labels:
        print("No clusters found. Try changing the parameters.")
        exit()


    print(f'Number of clusters: {num_clusters}\n')
    print(f'Best th = {best_th}, best size = {best_size}, best eps = {best_eps}\n')

    if args.out:
        for lab in unique_labels:
            if lab != -1:
                #cluster_points = [data[i] for i in range(len(data)) if labels[i] == lab]
                # Nou
                cluster_points = [tuple(data[i]) for i in range(len(data)) if labels[i] == lab]
                
                # Nou
                cluster_points = set(cluster_points)
                points_to_check = list(cluster_points)

                # Now, we will expand the cluster points to include neighbors within a radius of 3 residues
                #for point in cluster_points:
                #Nou
                for point in points_to_check:
                    for k in range(-3, 4):
                        for l in range(-3, 4):
                            if k**2 + l**2 <= 9: # Within a radius of 3 (residues)
                                new_point = [point[0] + k, point[1] + l]
                                #if new_point not in cluster_points:
                                #    cluster_points.append(new_point)
                                '''
                                if tuple(new_point) not in [tuple(cp) for cp in cluster_points]:
                                    cluster_points.append(new_point)
                                '''
                                if tuple(new_point) not in cluster_points and new_point[0] < len(zscore[0]) and new_point[0] >= 0 and new_point[1] < len(zscore[0]) and new_point[1] >= 0:
                                    cluster_points.add(tuple(new_point))
                string_points = [str(point[0]) + '-' + str(point[1]) for point in cluster_points]
                with open(args.out, 'a') as f:
                    f.write(' '.join(string_points) + '\n')


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
