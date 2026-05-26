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


def score(W, x, recompute_z=False):
    #x = str2int(x)

    # if length of sequence != length of model
    #if x.shape[-1] != len(mrf["v_idx"]):
    #    x = x[...,mrf["v_idx"]]
    
    # one hot encode
    #x = np.eye(states)[x] 
    
    # get non-gap positions
    #no_gap = 1.0 - x[...,-1]
    
    # remove gap from one-hot-encoding
    #x = x[...,:-1]
    
    # compute score
    #vw = mrf["v"] + np.tensordot(x,mrf["w"],2)
    #vw = np.tensordot(x,W,2)  # Ho fem sense la part de "v", perquè només volem mirar la contribució "coevolutiva". S'haruia de mirar si realment fa falta la part de "v", però crec que no.
    
    # ============================================================================================
    # Note, Z (the partition function) is a constant. In GREMLIN, V, W & Z are estimated using all
    # the original weighted input sequence(s). It is NOT recommended to recalculate Z with a 
    # different set of sequences. Given the common ERROR of recomputing Z, we include the option 
    # to do so, for comparison.
    # ============================================================================================
    #h = np.sum(np.multiply(x,vw),axis=-1)

    #####-------------- AIXÒ ESTÀ BÉ???? PROVAR SI ÉS EL MATEIX QUE AMB EL TENSORDOT ------------################
    not_gap = []
    #print(x)
    for k in range(len(x)):
        if x[k] != 20:
            not_gap.append(k)

    h = 0
    #print(not_gap)
    for i in not_gap:
        for j in not_gap:
            if i != j:
                h += W[i][j][x[i]][x[j]] 


    if recompute_z:
        z = logsumexp(vw, axis=-1)
        return np.sum((h-z), axis=-1)
    else:
        return np.sum(h, axis=-1)



def modify_seq(seq, region):
    # Com que al fer l'encode i passar de seq a int, i després traiem els gaps (lúltima columna), posarem gaps a totes les posicions que no estiguin a la regió, 
    # de forma que així només tindrem en compte les posicions que ens interessen (és a dir, traiem tots els residus que no estiguin a la regió d'interès).
    x = str2int(seq)
    for i in range(len(seq)):
        if i not in region:
            x[i] = 20
    return x



with open(args.clusters, 'r') as f:
    clusters = f.readlines()


# Read MRF
print("Abans de llegir el MRF\n")
lineas_w = []
longitud = 0
with open(args.mrf, 'r') as f:
    for line in f:
        #print(line)
        if line.startswith('V['):
            longitud += 1

        if line.startswith('W['):
            lineas_w.append(line)
print("Per aquí va bé\n")

#longitud = len(lineas_w)+1
w = [float('nan')] * longitud
for i in range(longitud):
    w[i] = [float('nan')] * longitud

print("Per aquí també va bé\n")

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

#w = np.array(w, dtype=float)
#print(w)

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
        if next_append == True:
            seqs.append(line_s)
            next_append = False
        else:
            seqs[-1] = seqs[-1] + line_s
del msa


for i, cluster in enumerate(clusters):
    cluster = cluster.strip()
    cluster = cluster.split()
    cluster = [int(x) for x in cluster]

    scores = []
    for seq in seqs:
        #print(cluster)
        x = modify_seq(seq, cluster)
        #print(x)
        scores.append(score(w, x))
    
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
    
############## També es podira mirar les seqüències que no són coevolutives i fer un fitxer amb aquestes seqüències, per treure contactes en comptes de posar-los ###############
