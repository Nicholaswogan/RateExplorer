import numpy as np
import utils
import fitting

def sat_CH3CN_nist(T):
    # NIST. Valid 288.3 - 362.3 K
    A,B,C = 4.27873, 1355.374, -37.853
    return utils.saturation_antoine(A, B, C, T)

def sat_CH3CN(T):
    # Yu et al. (2023)
    if T < 216.9:
        P = np.exp(21.542 - 5853.5/T-6.9785e-3*T-1.9407e-4*T**1.5)
    elif T >= 216.9 and T < 229.32:
        P = np.exp(23.405 - 5655/T - 0.72971*np.log(T) - 1.6902e-3*T - 1.9407e-4*T**1.5)
    else:
        P = 10.0**(4.6691 - 1583.4/(T - 15.263))
    
    return P*1e6

def sat_CH3CN_photochem(T):
    A, B, C = 3.31, -484.0, -397000.0
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'CH3CN'
    mu = 41.0519 # NIST
    T_triple = 229.32 # NIST
    T_critical = 545.0 # NIST
    T_ref = 250.0
    P_ref = sat_CH3CN(T_ref) # NIST

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    # Yu et al. (2023)
    T_data = np.linspace(100.0,545.0,20)
    P_data = np.array([sat_CH3CN(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_CH3CN_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_CH3CN_nist(T)/1e6 for T in entry['T']])
    entry['label'] = 'NIST'
    data.append(entry)

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()