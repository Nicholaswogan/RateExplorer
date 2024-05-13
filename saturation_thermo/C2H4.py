import numpy as np
import utils
import fitting

def sat_C2H4_nist(T):
    # NIST. Valid 149 to 188 K
    A,B,C = 3.87261, 584.146, -18.307
    return utils.saturation_antoine(A, B, C, T)

def sat_C2H4_Moses(T):
    # Moses 1992 Valid 77 - 155 K
    if T > 120:
        log10P = 6.74756 - 585/(T - 18.16)
        Psat = (10.0**log10P) # mm of Hg
    if T <= 120 and T > 104:
        log10P = 50.79 - 1703/T - 17.141*np.log10(T)
        Psat = (10.0**log10P) # mm of Hg
    if T <= 104 and T > 89:
        log10P = 8.724 -  901.6/(T - 2.555)
        Psat = (10.0**log10P) # mm of Hg
    if T <= 89:
        log10P = 1.5477 - 1038.1*(1/T - 0.011) + 16537*(1/T - 0.011)**2
        Psat = (10.0**log10P) # mm of Hg
        Psat = Psat/1e3
    Psat =  Psat/750.062 # mmHg * (bar/mmHg) = bar
    return Psat*1e6 # dynes/cm^2

def sat_C2H4_photochem(T):
    A, B, C = 8.41, -1130.0, -47100.0
    return np.exp(A + B/T + C/T**2.0) # bars

def main():
    gas_g = 'C2H4'
    mu = 28.0532 # NIST
    T_triple = 104.0 # NIST
    T_critical = 282.5 # NIST
    T_ref = 150.0
    P_ref = sat_C2H4_nist(T_ref) # NIST

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    # Moses et al 1992
    T_data = np.linspace(77,145,10)
    P_data = np.array([sat_C2H4_Moses(T) for T in T_data])

    # NIST
    T_data1 = np.linspace(149,188,10)
    P_data1 = np.array([sat_C2H4_nist(T) for T in T_data1])
    T_data = np.append(T_data,T_data1)
    P_data = np.append(P_data,P_data1)

    # NIST
    T_data = np.append(T_data,T_critical)
    P_data = np.append(P_data,50.6e6)

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_C2H4_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, 90, 1000)

if __name__ == '__main__':
    main()