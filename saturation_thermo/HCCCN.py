import numpy as np
import utils
import fitting

def sat_HCCCN(T):
    # Yu et al. (2023)
    if T < 280.0:
        P = 10.0**(9.0744 - 3046.5/(T + 36.176))
    else:
        P = 10.0**(4.6649 - 1470/T)
    return P*1e6

def sat_HCCCN_photochem(T):
    A, B, C = 1.68e+01, -5.09e+03, -9.79e-04
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'HCCCN'
    mu = 51.048
    T_triple = 280.0 # 
    T_critical = 527.0 # 
    T_ref = 300.0 # K
    P_ref = sat_HCCCN(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    T_data = np.linspace(160.0, 526.0, 20)
    P_data = np.array([sat_HCCCN(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_HCCCN_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()