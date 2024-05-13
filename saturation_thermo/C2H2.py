import numpy as np
import utils
import fitting

def sat_C2H2(T):
    # NIST. Valid 192.5 to 308 K
    if T < 214:
        A,B,C = 4.19598,699.53,-21.47
    else:
        A,B,C = 4.66141,909.079,7.947
    return utils.saturation_antoine(A, B, C, T)

def sat_C2H2_photochem(T):
    A, B, C = 10.7, -2070.0, 13300.0
    return np.exp(A + B/T + C/T**2.0) # bars

def sat_C2H2_Moses(T):
    # Moses 1992 valid 80 to 145 K
    log10P = 6.09748 - 1644.1/T + 7.42346*np.log10(1000/T)
    Psat = (10.0**log10P) # mm of Hg
    Psat =  Psat/750.062 # mmHg * (bar/mmHg) = bar
    return Psat*1e6 # dynes/cm^2

def main():
    gas_g = 'C2H2'
    mu = 26.0373
    T_triple = 192.4 # NIST
    T_critical = 308.3 # NIST
    T_ref = 200 # K
    P_ref = sat_C2H2(T_ref) # dynes/cm^2

    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)
    T_data = np.linspace(192.5,307,20)
    P_data = np.array([sat_C2H2(T) for T in T_data])

    # Handbook of chemistry and physics
    T_data = np.append([130,140,150,160,170,180,190],T_data)
    P_data = np.append(np.array([0.2,0.7, 2.6,7.8,20.6,49.0, 106.3])*(1e4/1),P_data)

    # Moses et al. (1992)
    T_data1 = np.array([[80,90,100,120]])
    P_data1 = np.array([sat_C2H2_Moses(T) for T in T_data1])
    T_data = np.append(T_data1,T_data)
    P_data = np.append(P_data1,P_data)

    guess = [4e9, 0, 6e9, 0]
    sol = h.optimize(guess, T_data, P_data)

    h.to_file('results/'+gas_g+'_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(T_data[0],T_critical,100)
    entry['P'] = np.array([sat_C2H2_photochem(T) for T in entry['T']])
    entry['label'] = 'Photochem_v0.5.4'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(80,145,100)
    entry['P'] = np.array([sat_C2H2_Moses(T)/1e6 for T in entry['T']])
    entry['label'] = 'Moses (1992)'
    data.append(entry)

    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, 120, 1000)

if __name__ == '__main__':
    main()