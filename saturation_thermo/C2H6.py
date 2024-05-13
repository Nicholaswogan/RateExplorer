import numpy as np
import utils
import fitting

def sat_C2H6_nist(T):
    # NIST. Valid 91 to 200 K
    if T < 135.7:
        A,B,C = 4.50706, 791.3, -6.422
    else:
        A,B,C = 3.93835,659.739,-16.719

    return utils.saturation_antoine(A, B, C, T)

def sat_C2H6_Moses(T):
    # Moses 1992 Valid 30 to 140 K
    if T < 90:
        log10P = 10.01 - 1085/(T - 0.561)
    else:
        log10P = 5.9366 - 1086.17/T + 3.83464*np.log10(1000/T)
    Psat = (10.0**log10P) # mm of Hg
    Psat =  Psat/750.062 # mmHg * (bar/mmHg) = bar
    return Psat*1e6 # dynes/cm^2

def sat_C2H6_photochem(T):
    A, B, C = 10.3, -1800.0, -14300.0
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'C2H6'
    mu = 30.0690 # NIST
    T_triple = 91.0 # NIST
    T_critical = 305.3 # NIST
    T_ref = 200.0
    P_ref = sat_C2H6_nist(T_ref) # NIST

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    # Moses et al 1992
    T_data = np.linspace(30,90,10)
    P_data = np.array([sat_C2H6_Moses(T) for T in T_data])

    # NIST
    T_data1 = np.linspace(100,200,10)
    P_data1 = np.array([sat_C2H6_nist(T) for T in T_data1])
    T_data = np.append(T_data,T_data1)
    P_data = np.append(P_data,P_data1)

    # NIST
    T_data = np.append(T_data,T_critical)
    P_data = np.append(P_data,49.0e6)

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_C2H6_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, T_triple-30, 1000)

if __name__ == '__main__':
    main()