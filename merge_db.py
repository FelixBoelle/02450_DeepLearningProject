# merge the CMR-low dimenstionality database with the OQMD structure files

from ase.db import connect
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from ase.io import write
from operator import itemgetter

def plot_n_atoms(lens):
    fig = sns.distplot(np.array(lens)[:,0])
    fig.set(xlabel='number of atoms')
    
# connect to databases
db_cmr = connect('mixdim.db')
db = connect('mixdim_DL.db')

entry_dic = {}
score_strings = []
dummy_line = db_cmr.get(id=1)
for key in vars(dummy_line)['_keys']:
    if 's_' in key:
        score_strings.append(key)
        entry_dic[key] = []

# visualize number of atoms first
lens = []
classes = []
for i in range(len(db_cmr)):
    i += 1
    row = db_cmr.get(id=i)
    
    if row.source == 'icsd':
        continue
    
    if i% 10000 == 0:
        print(i)
    # filter out all those weird overlapping guys, this is ICSD data with no values
    if row.volume == 1.0:
        continue
    n_atoms = len(row.positions)
    # filter out all unit cells with more than 200 atoms (less than 1%)
    if n_atoms > 200:
        continue
    
    score_value = 0
    # get dimenstionality class
    for score_string in score_strings:
        score_value_new = vars(row)[score_string]
        if score_value_new > score_value:
            dim_class = score_string
            score_value = score_value_new
        
    lens.append([len(db_cmr.get(id=i).positions),dim_class,i])
    entry_dic[dim_class].append([i,score_value])
    classes.append(dim_class)

# distribution over classes
for cl in score_strings:
    print(cl, ' number of entries: ',len(entry_dic[cl]))
    
# now pick 2000 random structures each
final_dic = {}
for key in vars(dummy_line)['_keys']:
    if 's_' in key:
        score_strings.append(key)
        final_dic[key] = sorted(entry_dic[key],key = itemgetter(1),reverse=True)[:2000]
        
# create final database, formula, positions, class
for clas in ['s_01','s_02','s_03','s_1','s_2','s_3']:
    for entry_id, score in final_dic[clas]:
        row = db_cmr.get(id=entry_id)
        db.write(row, formula = row.formula, score = score, dim_class = clas)
    