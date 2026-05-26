import matplotlib.pylab as plt
import argparse
from scipy import stats
import numpy as np

parser=argparse.ArgumentParser(description='Program description')
parser.add_argument('--matrix', required=True, help='The output file of gremlin_cpp')
parser.add_argument('--out', help='Output file (.txt) to write the zscore matrix')
parser.add_argument('--plot', help='If you do not want to plot it and just save the output, write an input here')
args=parser.parse_args()

matrix = args.matrix

def normalize(x):
  x = stats.boxcox(x - np.amin(x) + 1.0)[0]
  x_mean = np.mean(x)
  x_std = np.std(x)
  return((x-x_mean)/x_std)

def plot_mtx(apc, zscore):
  '''plot the mtx'''
  plt.figure(figsize=(15,5))
  for n, (key, data) in enumerate([("apc", apc), ("zscore", zscore)]):
    #plot
    plt.subplot(1,3,n+1)
    plt.title(key)
    if key == "zscore": 
       plt.imshow(data, cmap='Blues', vmin=1, vmax=3)
       plt.colorbar()
    else: 
       plt.imshow(data, cmap='Blues')
       plt.colorbar()
    plt.grid(False)
    
  plt.show()



with open(matrix, 'r') as f:
  lines = f.readlines()

columns = lines[-1].split()
num = int(columns[1]) + 1 # Posem +1 perquè comença per la posició 0
apc = np.zeros((num, num))
zscore_line = []
for line in lines[1:]:
  columns = line.split()
  apc[int(columns[0])][int(columns[1])] = float(columns[3])
  apc[int(columns[1])][int(columns[0])] = apc[int(columns[0])][int(columns[1])]
  zscore_line.append(float(columns[3]))

zscore_norm = normalize(zscore_line)
zscore = np.zeros((num, num))
idx = 0
for line in lines[1:]:
    columns = line.split()
    zscore[int(columns[0])][int(columns[1])] = zscore_norm[idx]
    zscore[int(columns[1])][int(columns[0])] = zscore[int(columns[0])][int(columns[1])]
    idx = idx + 1

if args.out:
    np.savetxt(args.out, zscore, delimiter=' ', fmt='%f')

if not args.plot:
  plot_mtx(apc, zscore)