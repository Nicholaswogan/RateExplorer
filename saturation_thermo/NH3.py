import numpy as np
import utils
import fitting

def sat_NH3(T):
    # Lodders & Fegley (1998), Table 1.20, Valid 160 K to 300 K
    if T < 195.4:
        a, b = 6.900, -1588
    else:
        a, b = 5.201, -1248
    return utils.saturation_lodders(a, b, T)

def sat_NH3_photochem(T):
    A = 14.8
    B = -5420.0
    C = 0.000711
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'NH3'
    mu = 17.031
    T_triple = 195.4 # Lodders & Fegley (1998), Table 1.20
    T_critical = 405.4 # NIST
    T_ref = 300 # K
    P_ref = sat_NH3(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)
    T_data = np.linspace(160,300,20)
    P_data = np.array([sat_NH3(T) for T in T_data])
    T_data = np.append(T_data, T_critical)
    P_data = np.append(P_data,113.0e6) # NIST gives P_critical

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(150,405,100)
    entry['P'] = np.array([sat_NH3_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem v0.5.4'
    data = [entry]

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, 160, 1000)

if __name__ == '__main__':
    main()