import pymol
from pymol import cmd
import argparse
import os
import json

parser=argparse.ArgumentParser(description='Program description')
parser.add_argument('--dir_1', required=True, help='Path of the PDB files directory (first)')
parser.add_argument('--dir_2', required=True, help='Path of the PDB files directory (second)')
parser.add_argument('--rmsd_file', required=True, help='Output file for the rmsd (.json)')
args=parser.parse_args()

dir_1 = args.dir_1
dir_2 = args.dir_2
rmsd_file = args.rmsd_file

pdbs_1 = [os.path.join(dir_1, nombre) for nombre in os.listdir(dir_1)]
pdbs_2 = [os.path.join(dir_2, nombre) for nombre in os.listdir(dir_2)]
pdbs = pdbs_1 + pdbs_2

rmsd_all_to_all={}
check = {pdb1: {pdb2: False for pdb2 in pdbs} for pdb1 in pdbs}
for pdb1 in pdbs:
    rmsd_all_to_all[pdb1]={}
    for pdb2 in pdbs:
        if check[pdb1][pdb2] == False:
            check[pdb1][pdb2] = True
            check[pdb2][pdb1] = True
            
            if pdb1 == pdb2:
                rmsd_all_to_all[pdb1][pdb2] = 0
            else:
            # Open the PDB files
                cmd.load(f'{pdb1}', 'mol1')
                cmd.load(f'{pdb2}', 'mol2')

                #If there is any HETATM that could give us further problems when performing NMA (usually an "unatached" CA), they will be removed
                cmd.select("nhmol1", '%s and not %s'%("mol1", "hetatm"))
                cmd.select("Ca1", '%s and name %s'%("nhmol1", "CA"))
                cmd.select("nhmol2", '%s and not %s'%("mol2", "hetatm"))
                cmd.select("Ca2", '%s and name %s'%("nhmol2", "CA"))
                
                #print("\n\n", pdb1, " ", pdb2, "\n\n")
                rmsd_all_to_all[pdb1][pdb2] = cmd.align("Ca1", "Ca2", cycles=0)[0]
                cmd.delete("all")
                cmd.reinitialize()
        else:
            rmsd_all_to_all[pdb1][pdb2] = rmsd_all_to_all[pdb2][pdb1]



with open(rmsd_file, "w") as archivo:
    json.dump(rmsd_all_to_all, archivo)
