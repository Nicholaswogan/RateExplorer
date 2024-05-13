import numpy as np
import utils
import fitting

def sat_CH4(T):
    # Yu et al. (2023)
    if T < 90.686:
        P = 10.0**(4.31972- 451.64/(T - 4.66))
    else:
        P = 10.0**(3.96867 - 435.453/(T - 1.789))
    return P*1e6

def main():
    gas_g = 'CH4'
    mu = 16.04
    T_triple = 90.686 # Yu et al. (2023)
    T_critical = 190.564 # Yu et al. (2023)
    T_ref = 100.0 # K
    P_ref = sat_CH4(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    # Yu et al. (2023)
    T_data = np.linspace(30.0, 190.0, 20)
    P_data = np.array([sat_CH4(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    fitting.plot_sat(h, gas_g, T_data, P_data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()