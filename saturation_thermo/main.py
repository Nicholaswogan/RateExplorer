import numpy as np
from matplotlib import pyplot as plt
import utils
from photochem.utils._format import yaml, flowmap, blockseqtrue, MyDumper

SPECIES = [
    'H2O', 'CO2', 'S8', 'H2SO4', # Saturates in rocky planet atmospheres
    'NH3', 'N2O', # Saturates on Jupiter (not N2O?)
    'C2H2', 'C2H4', 'C2H6', 'CH3CN', 'HCCCN', 'HCN', 'CH4' # Saturates on Titan
]

def plot_all():
    
    plt.rcParams.update({'font.size': 12})
    fig,ax = plt.subplots(1,1,figsize=[6,5])

    for i,sp in enumerate(SPECIES):
        h = utils.SaturationProperties()
        h.init2('results/'+sp+'_sat.yaml')
        TT = np.linspace(50,h.T_critical,100)
        val = np.array([h.sat_pressure(T) for T in TT])/1e6

        if i > 9:
            ax.plot(TT, val, label=sp, ls='--')
        else:
            ax.plot(TT, val, label=sp)

    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')
    ax.legend(ncol=1,bbox_to_anchor=(1.03, 1.0), loc='upper left')
    ax.grid(alpha=0.4)
    ax.set_ylim(1e-20,1e3)
    plt.savefig('figures/all.pdf',bbox_inches='tight')
    plt.close()

def main():
    for sp in SPECIES:
        exec('import '+sp)
        exec(sp+'.main()')

def collect_files():
    out1 = []
    tmp1 = []
    for sp in SPECIES:
        with open('results/'+sp+'_sat.yaml') as f:
            sat = yaml.load(f,Loader=yaml.Loader)
        with open('results/'+sp+'_thermo.yaml') as f:
            thermo = yaml.load(f,Loader=yaml.Loader)[0]
        out = {}
        out['name'] = sp+'aer'
        out['composition'] = flowmap(thermo['composition'])
        out['density'] = 1.0
        out['optical-properties'] = "none"
        out['formation'] = 'saturation'
        out['gas-phase'] = sp
        sat['parameters'] = flowmap(sat['parameters'])
        sat['vaporization'] = flowmap(sat['vaporization'])
        sat['sublimation'] = flowmap(sat['sublimation'])
        sat['super-critical'] = flowmap(sat['super-critical'])
        out['saturation'] = sat

        out1.append(out)

        tmp = {}
        tmp['name'] = sp+'aer'
        tmp['composition'] = flowmap(thermo['composition'])
        tmp['condensate'] = True
        thermo['thermo']['temperature-ranges'] = blockseqtrue(thermo['thermo']['temperature-ranges'])
        thermo['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(thermo['thermo']['data'])]
        tmp['thermo'] = thermo['thermo']
        tmp1.append(tmp)
    
    with open('results/condensates.yaml','w') as f:
        yaml.dump(out1,f,Dumper=utils.CustomDumper,sort_keys=False,width=70)

    with open('results/condensates_shomate.yaml','w') as f:
        yaml.dump(tmp1,f,Dumper=utils.CustomDumper,sort_keys=False,width=70)     
            
if __name__ == '__main__':
    main()
    plot_all()
    collect_files()