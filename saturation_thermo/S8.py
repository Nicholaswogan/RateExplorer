import numpy as np
import utils
import fitting

def sat_sulfur_vap(a1, a2, T):
    P = 10.0**(a1 - a2/T)
    P = P*1.01325
    P = P*1e6
    return P

def sat_S2_vap(T):
    a1, a2 = 7.0240, 6091.2
    return sat_sulfur_vap(a1, a2, T)

def sat_S3_vap(T):
    a1, a2 = 6.3428, 6202.2
    return sat_sulfur_vap(a1, a2, T)

def sat_S4_vap(T):
    a1, a2 = 6.0028, 6047.5
    return sat_sulfur_vap(a1, a2, T)

def sat_S5_vap(T):
    a1, a2 = 5.1609, 4714.8
    return sat_sulfur_vap(a1, a2, T)

def sat_S6_vap(T):
    a1, a2 = 4.8039, 3814.1
    return sat_sulfur_vap(a1, a2, T)

def sat_S7_vap(T):
    a1, a2 = 5.2127, 4113.6
    return sat_sulfur_vap(a1, a2, T)

def sat_S8_vap(T):
    a1, a2 = 4.1879, 3269.1
    return sat_sulfur_vap(a1, a2, T)

def sat_all_S_vap(T):
    P = sat_S2_vap(T) + sat_S3_vap(T) + sat_S4_vap(T) + sat_S5_vap(T) + sat_S6_vap(T)
    + sat_S7_vap(T) + sat_S8_vap(T)
    return P

def sat_S8(T):
    # Lyons (2008)
    if T > 392.0:
        P = sat_S8_vap(T)
    else:
        P_i_liq = sat_S8_vap(T)

        a1, a2 = 8.4832, 5082.0
        P_sol = sat_sulfur_vap(a1, a2, T)

        P_liq = sat_all_S_vap(T)

        P = P_i_liq*(P_sol/P_liq)/10

    return P

def sat_S8_photochem(T):
    A, B, C = 8.09, -4941.34, -940367.8
    return np.exp(A + B/T + C/T**2.0) # bars

def sat_S8_rimmer(T):
    return sat_S8_vap(T)

def sat_S8_zahnle(T):
    if T < 413:
        P = np.exp(20 - 11800/T)
    else:
        P = np.exp(9.6 - 7510/T)
    return P*1e6

def main():
    gas_g = 'S8'
    mu = 256.520
    T_triple = 392.0 # Lyons (2008)
    T_critical = 675.0 # Guess
    T_ref =  400.0
    P_ref = sat_S8(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    T_data = np.linspace(220,675,20)
    P_data = np.array([sat_S8(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S8_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem v0.5.4'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S8_rimmer(T)/1e6 for T in entry['T']])
    entry['label'] = 'Rimmer (2021)'
    data.append(entry)

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_S8_zahnle(T)/1e6 for T in entry['T']])
    entry['label'] = 'Zahnle (2016)'
    data.append(entry)

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()