import os
from utils import yaml
import utils

def runall():
    species = [a for a in os.listdir('.') if '.py' in a and a not in ['main.py','template.py','utils.py']]
    for sp in species:
        sp1 = sp.strip('.py')
        print(sp1)
        exec('import '+sp1)
        exec(sp1+'.main()')
        print()

def collect_citations():
    species = list(set([a.replace('.h5','').replace('.yaml','') for a in os.listdir('results') if a != '.gitignore' and a != 'metadata.yaml']))

    citations = {}
    for i,sp in enumerate(species):
        with open('results/'+sp+'.yaml','r') as f:
            dat = yaml.load(f,Loader=yaml.Loader)
        dat = utils.format_citation(sp,dat)
        citations[sp] = dat[sp]
    
    with open('results/metadata.yaml','w') as f:
        yaml.dump(citations, f, utils.MyDumper, sort_keys=False)

def main():
    runall()
    collect_citations()

if __name__ == '__main__':
    main()

    


