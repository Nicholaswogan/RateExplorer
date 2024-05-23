import numpy as np
import utils
import fitting

def sat_H2SO4_rimmer(T):
    log10P = 4.4753 - 3229.0/(T + 7.1192e5) - 3.1723e6/T**2 + 4.0832e8/T**3 - 2.0312e10/T**4
    return 10.0**log10P

def sat_H2SO4(T):
    # Dai et al. (2022)
    Tc = 905
    T0 = 360
    lnP = 16.259 - 10156/T0 + 10156*(-1/T + 1/T0 + (0.38/(Tc-T0))*(1 + np.log(T0/T) - T0/T))
    return np.exp(lnP)*1e6

def main():
    gas_g = 'H2SO4'
    mu = 98.079
    T_triple = 360.0 # Dai et al. (2022)
    T_critical = 905.0 # Dai et al. (2022)
    T_ref = 400.0 # K
    P_ref = sat_H2SO4(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)
    h.dlog10P_critical = 1
    
    # 
    T_data = np.linspace(200.0, 900.0, 20)
    P_data = np.array([sat_H2SO4(T) for T in T_data])

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    fitting.plot_sat(h, gas_g, T_data, P_data)
    fitting.fit_thermo(gas_g, T_triple-30, 1500)

if __name__ == '__main__':
    main()