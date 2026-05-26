import numpy as np
import os
from scipy.special import logsumexp
from scipy.stats import zscore
import argparse

parser = argparse.ArgumentParser(description='Calculate the coevolutionary score of all sequences using a GREMLIN model, and cluster them using by "high coevolution".')
parser.add_argument('--mrf', required=True, help='Path to the GREMLIN model file (MRF).')
parser.add_argument('--msa', required=True, help='Path to the sequences file (fasta).')
parser.add_argument('--clusters', required=True, help='Path to the clusters file (.txt).')
parser.add_argument('--code', required=True, help='Code of the protein to which the MSA corresponds to (e.g., 1RGS)')
parser.add_argument('--out_dir', required=True, help='Output directory for the MSA clusters.')
parser.add_argument('--th', required=True, help='Threshold for z-score to consider a residue as coevolving.')
args = parser.parse_args()

th = float(args.th)


###############
## FUNCTIONS
###############

alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)

# map amino acids to integers (A->0, R->1, etc)
a2n = dict((a,n) for n,a in enumerate(alphabet))
aa2int = lambda x: a2n.get(x,a2n['-'])


def str2int(x):
    '''convert a list of strings into list of integers'''
    # Example: ["ACD","EFG"] -> [[0,4,3], [6,13,7]]
    x = list(x)
    x = np.array(x)
    '''
    if x.dtype.type is np.str_:
        if x.ndim == 0: return np.array([aa2int(aa) for aa in x])
        else: return np.array([[aa2int(aa) for aa in seq] for seq in x])
    else:
        return x
    '''
    return np.array([aa2int(aa) for aa in x])


def score(W, x, pairs):
    h = 0
    for pair in pairs:
        i = pair[0]
        j = pair[1]
        if x[i] != 20 and x[j] != 20: # We only compute the score for non-gap positions
            h += W[i][j][x[i]][x[j]]
    return h


if __name__ == "__main__":
    with open(args.clusters, 'r') as f:
        clusters = f.readlines()

    # Read MRF
    lineas_w = []
    longitud = 0
    with open(args.mrf, 'r') as f:
        for line in f:
            #print(line)
            if line.startswith('V['):
                longitud += 1

            if line.startswith('W['):
                lineas_w.append(line)

    w = [float('nan')] * longitud
    for i in range(longitud):
        w[i] = [float('nan')] * longitud


    for line in lineas_w:
        # Eliminar saltos de línea y espacios extra
        line = line.strip()
        # Dividir la línea en partes (por espacios)
        parts = line.split()
        
        # Ignorar el primer elemento que es 'W[i][j]' y tomar solo los números
        numbers = parts[1:]  # Tomar solo los números a partir de la segunda parte
        numbers_float = np.array(numbers, dtype=float)
        matriz = np.reshape(numbers_float, (21, 21))

        indices = parts[0]
        fila, columna = map(int, indices.strip("W[]").split("]["))
        w[fila][columna] = matriz
        w[columna][fila] = matriz
    del lineas_w


    with open(args.msa, 'r') as f:
        msa = f.readlines()

    # El següent for és perquè el MSA no es una seqüència per línea (perquè hem utilitzat mafft), i per tant hem de veure com sel·leccionem les seqüències, que serà una seqüència nova cada cop que hi hagi un ">"
    next_append = False # Ens servirà per saber si la línea anterior era una referència, per saber si hem de fer append o +
    refs = []
    seqs = []
    for line in msa:
        line_s = line.strip()
        if line_s.startswith('>'):
            refs.append(line_s)
            next_append = True
        else:
            if next_append:
                seqs.append(line_s)
                next_append = False
            else:
                seqs[-1] = seqs[-1] + line_s
    del msa


    for i, cluster in enumerate(clusters):
        cluster = cluster.strip()
        cluster = cluster.split(' ')
        points = []
        for c in cluster:
            c = c.split('-')
            if len(c) != 2:
                print(f'Error in cluster {i} at point {c}')
                exit(1)
            points.append((int(c[0]), int(c[1])))

        scores = []
        for seq in seqs:
            x = str2int(seq)
            scores.append(score(w, x, points))
        
        zscores = zscore(scores)
        indices = []
        for j in range(len(zscores)):
            if zscores[j] >= th:
                indices.append(j)
        
        # Write the sequences that are above the threshold
        out_name = f'OUT_{args.code}_{i:03}.fasta'
        out = os.path.join(args.out_dir, out_name)
        with open(out, 'a') as o:
            o.write(f'{refs[0]}\n')
            o.write(f'{seqs[0]}\n')
            for j in indices:
                if j != 0:
                    o.write(f'{refs[j]}\n')
                    o.write(f'{seqs[j]}\n')