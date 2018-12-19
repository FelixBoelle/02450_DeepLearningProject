# after preparing the databases, now useful fingerprints should be identified and extracted

from ase.db import connect
from ase.data import covalent_radii
import numpy as np
import matplotlib.pyplot as plt
from ase.visualize import view
from ase.build import bulk
from operator import itemgetter

def get_all_k_values(atoms):
    ''' calculate the k-value according to the covalent bond with its 20 nearest
    neighbors '''
    
    ks = []
    for i, atom_i in enumerate(atoms):
        ks_temp = []
        for j, atom_j in enumerate(atoms):
            
            if i == j: # filter out 0 distance to itself
                continue

            d_ij = distances[i][j]
            k = get_k_value(atom_i.number, atom_j.number, d_ij)     
            ks_temp.append(np.around(k, decimals = 3))
        
        # filter stuff that has overlapping atoms
        if any(i > 1.3 for i in ks_temp):
            print("Some elements overlap in here id: {}".format(atoms.get_chemical_formula()))
            return "overlap"

        ks.append(sorted(ks_temp, reverse = True)[:n_nearest_neighbors])
    
    return ks

def get_k_value(n_i,n_j,d_ij):
    ''' calculate k-value given by low dimensionality paper
    
    input:
        n_i: atomic number atom i
        n_j: atomix numer atom j
        d_ij: distance in real space atom i-j in angstrom
    returns:
        k-value '''
        
    k = (covalent_radii[n_i] + covalent_radii[n_j])/d_ij
    
    return k

def padding_no_order(ks, n_max):
    ''' Padding applied by just repeating pattern and using n_max as cut-off '''
    
    ks += ks
    ks = ks[:n_max]
    return ks

#
def padding_with_order(ks, n_max):
    ''' Padding applied by just repeating pattern and using n_max as cut-off '''
    
    ks += ks
    ks = ks[:n_max]
    mean_values = sorted([[np.mean(k),k] for i,k in enumerate(ks)], 
                           key = itemgetter(0))

    ks = [val[1] for val in mean_values]
        
    return ks

db = connect('mixdim_DL.db')
db1 = connect('s_1.db')
db2 = connect('s_2.db')
db3 = connect('s_3.db')
db4 = connect('s_01.db')
db5 = connect('s_02.db')
db6 = connect('s_03.db')

# how many nearest neighbors to include
n_max = 200
n_nearest_neighbors = 200
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

# the different classes used for the dimensionality classifciation
classes = ['s_1','s_2','s_3','s_01','s_02','s_03']

for n_cla,cla in enumerate(classes):
    print('_________{}___________'.format(cla))
    for row in db.select(dim_class = cla):
        atoms = row.toatoms()
        print(atoms.get_chemical_formula())
        
        if atoms.get_chemical_formula() == 'NiO2':
            continue
       
        # repeat cell in a fair way to get at least 200 atoms in unit cell
        n_repeat = int(n_max/len(atoms))
        padding = n_max - n_repeat*len(atoms)
        for i in range(n_max):
           
            if len(atoms) * 2 > n_max:
                break
            
            cell = atoms.get_cell()                                                                                                                                                                                                                                 
            x = np.linalg.norm(atoms.get_cell()[0,:])
            y = np.linalg.norm(atoms.get_cell()[1,:])                                                                                                                                                                                                                                                           
            z = np.linalg.norm(atoms.get_cell()[2,:])
        
            if x <= z and x <= y:
                atoms = atoms.repeat((2,1,1))
                continue
            if y <= z and y <= x:
                atoms = atoms.repeat((1,2,1))
                continue
            if z <= y and z <= x:
                atoms = atoms.repeat((1,1,2))
                continue
        
        distances = atoms.get_all_distances(mic = True)
        
        # all k-values in structure first
        ks = get_all_k_values(atoms)
        if ks == "overlap":
            print('overlap for {} {}'.format(row.id, cla))
            continue
        
        # padding by repeating pattern
        ks = padding_no_order(ks, n_max)
        #ks = padding_with_order(ks, n_max)
        # visualize it using some sort of heat map, black and white is fine
        if 1:
            plt.figure(figsize = (20,8))
            
            plt.imshow(np.array(ks),cmap = 'inferno')
            view(atoms)
        
        # finally all entries are saved in separate databases, one for each class
        # this should provide greater flexibility
        if n_cla == 0: db2.write(atoms,data = {'ks':ks});
        if n_cla == 1: db3.write(atoms,data = {'ks':ks});
        if n_cla == 2: db3.write(atoms,data = {'ks':ks});
        if n_cla == 3: db4.write(atoms,data = {'ks':ks});
        if n_cla == 4: db5.write(atoms,data = {'ks':ks});
        if n_cla == 5: db6.write(atoms,data = {'ks':ks});