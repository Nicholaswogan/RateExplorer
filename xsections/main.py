import os
import numpy as np
import h5py
from utils import yaml
import utils

def runall():
    species = [a for a in os.listdir('.') if '.py' in a and a not in ['main.py','template.py','utils.py']]
    species.sort()

    for sp in species:
        sp1 = sp.strip('.py')
        print(sp1)
        exec('import '+sp1)
        exec(sp1+'.main()')
        print()

def collect_citations():
    species = list(set([a.replace('.h5','') for a in os.listdir('results') if '.h5' in a]))
    species.sort()

    citations = {}
    for i,sp in enumerate(species):
        with open('results/'+sp+'.yaml','r') as f:
            dat = yaml.load(f,Loader=yaml.Loader)
        dat = utils.format_citation(sp,dat)
        citations[sp] = dat[sp]
    
    with open('results/metadata.yaml','w') as f:
        yaml.dump(citations, f, utils.MyDumper, sort_keys=False)

    for i,sp in enumerate(species):
        os.remove('results/'+sp+'.yaml')

def make_bins_file():
    wavl = []

    with open('Zahnle_Kevin_data/SW_neutrals_for_stays.DAT','r') as f:
        lines = f.readlines()

    for line in lines[2:]:
        wv = float(line.split()[1])/10
        if np.isclose(wv,175.4):
            break
        wavl.append(wv)
        
    with open('Zahnle_Kevin_data/photo_new_new.dat','r') as f:
        lines = f.readlines()
    for line in lines[2:]:
        if len(line) < 5:
            break
        tmp = line.split()[1]
        wavl.append(float(tmp.split('-')[0])/10)
    wavl.append(8550.0/10)

    wavl = wavl[5:]

    while True:
        if wavl[-1]+10.0 > 1100.0:
            break
        wavl.append(wavl[-1]+10.0)
    wavl[-1] = 1100.0

    wavl = np.array(wavl)

    with h5py.File('results/bins.h5','w') as f:
        dset = f.create_dataset("wavl", wavl.shape, 'f')
        dset[:] = wavl

def main():
    runall()
    collect_citations()
    make_bins_file()

if __name__ == '__main__':
    main()

    


