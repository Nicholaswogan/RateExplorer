import numpy as np
import utils
import fitting
import S8

def sat_S2(T):
    # Lyons (2008)
    if T > 392.0:
        P = S8.sat_S2_vap(T)
    else:
        P_i_liq = S8.sat_S2_vap(T)

        a1, a2 = 8.4832, 5082.0
        P_sol = S8.sat_sulfur_vap(a1, a2, T)

        P_liq = S8.sat_all_S_vap(T)

        P = P_i_liq*(P_sol/P_liq)/10

    return P

def sat_S2_photochem(T):
    A, B, C = 14.5, -11318.6, -984857.9
    return np.exp(A + B/T + C/T**2.0) # bars

def sat_S2_rimmer(T):
    return S8.sat_S2_vap(T)

def sat_S2_zahnle(T):
    if T < 413:
        P = np.exp(27 - 18500/T)
    else:
        P = np.exp(16.1 - 14000/T) 
    return P*1e6

def main():
    gas_g = 'S2'
    mu = 64.13
    T_triple = 392.0 # Lyons (2008)
    T_critical = 675.0 # Guess
    T_ref = 400.0
    P_ref = sat_S2(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    T_data = np.linspace(220,675,20)
    P_data = np.array([sat_S2(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S2_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem v0.5.4'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S2_rimmer(T)/1e6 for T in entry['T']])
    entry['label'] = 'Rimmer (2021)'
    data.append(entry)

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S2_zahnle(T)/1e6 for T in entry['T']])
    entry['label'] = 'Zahnle (2016)'
    data.append(entry)

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()