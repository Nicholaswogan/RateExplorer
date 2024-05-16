import numpy as np
from matplotlib import pyplot as plt
import utils

SPECIES = ['H2O','CO2','NH3','C2H2','C2H4','C2H6','CH3CN',
           'N2O','HCCCN','HCN','CH4','S2','S3','S4','S8','H2SO4']

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
    
if __name__ == '__main__':
    main()
    plot_all()