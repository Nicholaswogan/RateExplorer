import numpy as np
import utils
import fitting

def sat_HCN(T):
    # Lodders & Fegley (1998), Table 1.20
    if T < 259:
        a, b = 6.747, -1944.5
    else:
        a, b = 5.037, -1504.5
    return utils.saturation_lodders(a, b, T)

def main():
    gas_g = 'HCN'
    mu = 27.0253
    T_triple = 259.0 # Lodders & Fegley (1998), Table 1.20
    T_critical = 456.85 # Matweb
    T_ref = 290.0 # K
    P_ref = sat_HCN(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)
    T_data = np.linspace(100,299,20)
    P_data = np.array([sat_HCN(T) for T in T_data])
    T_data = np.append(T_data, T_critical)
    P_data = np.append(P_data,53.9e6) # Matweb

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    fitting.plot_sat(h, gas_g, T_data, P_data)
    fitting.fit_thermo(gas_g, 140, 1000)

if __name__ == '__main__':
    main()