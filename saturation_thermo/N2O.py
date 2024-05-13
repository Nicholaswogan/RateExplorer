import numpy as np
import utils
import fitting

def sat_N2O_nist(T):
    # NIST. Valid 129.8 - 187.7 K
    A,B,C = 4.37799, 621.077, -44.659
    return utils.saturation_antoine(A, B, C, T)

def sat_N2O_vap(T):
    a1, a2, a3, a4 = -6.8657, 1.9373, -2.6440, 0.0387
    tau = (1 - T/309.56)
    P = 72.38e6*np.exp((a1*tau + a2*tau**1.5 + a3*tau**2.5 + a4 *tau**5)/(T/309.56))
    return P

def sat_N2O(T):
    # Ferreira and Lobo (2009)

    if T > 182.33:
        P = sat_N2O_vap(T)
    else:
        Tt = 182.33
        e1, e2 = -6.6551, -9.8076
        f = 1.0364
        theta = (1.0 - T/Tt)
        Pt = sat_N2O_vap(Tt)/1e4 # dynes/cm^2 * (kPa / 10000 dynes/cm^2) = kPa
        P = np.exp(np.log(Pt) + (e1*theta + e2*theta**f)/(T/Tt)) # kPa
        P = P*1e4 # to dynes/cm^2
        
    return P

def sat_N2O_photochem(T):
    A, B, C = 1.62e+01, -2.97e+03, -7.01e-04
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'N2O'
    mu = 44.0128 # NIST
    T_triple = 182.33 # Ferreira and Lobo (2009)
    T_critical = 309.56 # NIST
    T_ref = 185.0
    P_ref = sat_N2O_nist(T_ref) # NIST

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    # Ferreira and Lobo (2009)
    T_data = np.linspace(70.0,309.0,20)
    P_data = np.array([sat_N2O(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_N2O_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(129.8, 187.7, 100)
    entry['P'] = np.array([sat_N2O_nist(T)/1e6 for T in entry['T']])
    entry['label'] = 'NIST'
    data.append(entry)

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()