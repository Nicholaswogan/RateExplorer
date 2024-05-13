
import numpy as np
from matplotlib import pyplot as plt
import utils
from scipy import optimize
import copy
from photochem.utils._format import yaml

def plot_sat(h, gas_g, T_data, P_data, data = None):
    plt.rcParams.update({'font.size': 12})
    fig,axs = plt.subplots(1,2,figsize=[10,4])

    ax = axs[0]
    TT = np.linspace(T_data[0],h.T_critical+1,1000)
    P = np.array([h.sat_pressure(h.A_v,h.B_v,h.A_s,h.B_s,h.A_c,h.B_c,T) for T in TT])/1e6
    ax.plot(T_data, P_data/1e6,marker='o',ms=3,ls='',label='Data')
    ax.plot(TT, P,'--',label='Fit')
    if data is not None:
        for d in data:
            ax.plot(d['T'], d['P'],':',label=d['label'])
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')
    ax.axvline(h.T_triple,c='k',ls=':',lw=0.7)
    ax.axvline(h.T_critical,c='k',ls=':',lw=0.7)
        
    ax = axs[1]
    TT = np.linspace(100,h.T_critical,100)
    val = np.array([h.latent_heat(h.A_v,h.B_v,h.A_s,h.B_s,h.A_c,h.B_c,T) for T in TT])
    ax.plot(TT, val,'--',label='Fit')
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Latent heat (erg/g)')
    ax.set_xlabel('Temperature (K)')
    ax.axvline(h.T_triple,c='k',ls=':',lw=0.7)
    ax.axvline(h.T_critical,c='k',ls=':',lw=0.7)

    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    plt.savefig('figures/'+gas_g+'_sat.pdf',bbox_inches='tight')
    plt.close()

def Psat_coeffs(r, gas_g, coeffs_c, T):
    Gs = utils.gibbs_energy_shomate(coeffs_c, T)
    Gg = r.gibbs_energy(gas_g,T)
    return utils.Psat_eqn(Gs, Gg, T)

def objective_(x, r, gas_g, Psat_data, T_data):
    coeffs_c = np.zeros(7)
    coeffs_c[0] = x[0]
    coeffs_c[5] = x[1]
    coeffs_c[6] = x[2]
    Psat = np.array([Psat_coeffs(r, gas_g, coeffs_c, T) for T in T_data])
    fval1 = np.log10(Psat) - np.log10(Psat_data)
    return fval1

def fit_thermo(gas_g, T_low, T_high):
    r = utils.ReactionExplorer('thermodata121.yaml')

    sat_fcn = utils.SaturationProperties()
    sat_fcn.init2('results/'+gas_g+'_sat.yaml')

    T_data = np.linspace(T_low,sat_fcn.T_triple,100)
    Psat_data = np.array([sat_fcn.sat_pressure(T) for T in T_data])/1e6 # bars
    x_init = [
        r.heat_capacity(gas_g,np.mean(T_data)),
        r.enthalpy(gas_g,np.mean(T_data))/1e3,
        r.entropy(gas_g,np.mean(T_data))
    ]
    sol1 = optimize.root(objective_, x_init, method='lm', args=(r, gas_g, Psat_data, T_data))
    if not sol1.success:
        raise Exception()

    T_data = np.linspace(sat_fcn.T_triple,sat_fcn.T_critical,100)
    Psat_data = np.array([sat_fcn.sat_pressure(T) for T in T_data])/1e6 # bars
    x_init = [
        r.heat_capacity(gas_g,np.mean(T_data)),
        r.enthalpy(gas_g,np.mean(T_data))/1e3,
        r.entropy(gas_g,np.mean(T_data))
    ]
    sol2 = optimize.root(objective_, x_init, method='lm', args=(r, gas_g, Psat_data, T_data))
    if not sol2.success:
        raise Exception()

    T_data = np.linspace(sat_fcn.T_critical,T_high,100)
    Psat_data = np.array([sat_fcn.sat_pressure(T) for T in T_data])/1e6 # bars
    x_init = [
        r.heat_capacity(gas_g,np.mean(T_data)),
        r.enthalpy(gas_g,np.mean(T_data))/1e3,
        r.entropy(gas_g,np.mean(T_data))
    ]
    sol3 = optimize.root(objective_, x_init, method='lm', args=(r, gas_g, Psat_data, T_data))
    if not sol3.success:
        raise Exception()

    thermo = {}
    thermo['model'] = 'Shomate'
    thermo['temperature-ranges'] = [0.0,sat_fcn.T_triple,sat_fcn.T_critical,6000]
    thermo['data'] = []
    thermo['data'].append([float('%e'%sol1.x[0]),0.0,0.0,0.0,0.0,float('%e'%sol1.x[1]),float('%e'%sol1.x[2])])
    thermo['data'].append([float('%e'%sol2.x[0]),0.0,0.0,0.0,0.0,float('%e'%sol2.x[1]),float('%e'%sol2.x[2])])
    thermo['data'].append([float('%e'%sol3.x[0]),0.0,0.0,0.0,0.0,float('%e'%sol3.x[1]),float('%e'%sol3.x[2])])
    entry = copy.deepcopy(r.mech['species'][gas_g])
    entry['name'] = gas_g+'[c]'
    entry['condensate'] = True
    entry['triple-temperature'] = sat_fcn.T_triple
    entry['critical-temperature'] = sat_fcn.T_critical
    entry['thermo'] = thermo
    del entry['note']

    with open('results/'+gas_g+'_thermo.yaml','w') as f:
        yaml.dump([utils.Format_species(entry)],f,Dumper=yaml.Dumper,sort_keys=False,width=70)

    # plot

    plt.rcParams.update({'font.size': 12})
    fig,axs = plt.subplots(2,4,figsize=[20,8])

    ax = axs[0,0]
    TT = np.linspace(50,sat_fcn.T_triple,100)
    val = np.array([sat_fcn.sat_pressure(T) for T in TT])/1e6
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.Psat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')

    ax = axs[0,1]
    TT = np.linspace(sat_fcn.T_triple-5,sat_fcn.T_triple+5,100)
    val = np.array([sat_fcn.sat_pressure(T) for T in TT])/1e6
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.Psat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')
    ax.axvline(sat_fcn.T_triple,c='k',ls=':',lw=0.7)

    ax = axs[0,2]
    TT = np.linspace(sat_fcn.T_triple,sat_fcn.T_critical,100)
    val = np.array([sat_fcn.sat_pressure(T) for T in TT])/1e6
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.Psat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')

    ax = axs[0,3]
    TT = np.linspace(sat_fcn.T_critical-10,6000,100)
    val = np.array([sat_fcn.sat_pressure(T) for T in TT])/1e6
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.Psat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.set_yscale('log')

    ax = axs[1,0]
    TT = np.linspace(100,sat_fcn.T_critical,100)
    val = np.array([sat_fcn.latent_heat(T)*(sat_fcn.mu/1)*(1/1e7) for T in TT])
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.latent_heat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Latent heat (J/mol)')
    ax.set_xlabel('Temperature (K)')
    
    ax = axs[1,1]
    TT = np.linspace(sat_fcn.T_critical-10,6000,100)
    val = np.array([sat_fcn.latent_heat(T)*(sat_fcn.mu/1)*(1/1e7) for T in TT])
    ax.plot(TT, val,'--',label='Truth')
    val = np.array([utils.latent_heat_eval(thermo, r.mech['species'][gas_g]['thermo'], T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Latent heat (J/mol)')
    ax.set_xlabel('Temperature (K)')

    ax = axs[1,2]
    TT = np.linspace(60,sat_fcn.T_critical,100)
    val = np.array([utils.heat_capacity_eval(thermo, T) for T in TT])
    ax.plot(TT, val,'-',label='Fit', lw=1)
    ax.legend()
    ax.grid(alpha=0.4)
    ax.set_ylabel('Heat capacity (J/(mol K))')
    ax.set_xlabel('Temperature (K)')

    ax = axs[1,3]
    ax.set_visible(False)

    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    plt.savefig('figures/'+gas_g+'_thermo.pdf',bbox_inches='tight')
    plt.close()
