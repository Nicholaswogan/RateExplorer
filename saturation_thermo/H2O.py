import numpy as np
from matplotlib import pyplot as plt
import utils
from scipy import optimize
import copy
from photochem.utils._format import yaml
import fitting

def get_data():

    T, PH2O, a, a = np.loadtxt('data/H2O_highT.txt').T
    T = T + 273.15
    PH2O = PH2O/1e2*1e6

    T1, PH2O1 = np.loadtxt('data/H2O_lowT.txt')
    T1 = T1 + 273.15
    PH2O1 = PH2O1*10

    T_SOL = np.append(T1,T)
    PH2O_SOL = np.append(PH2O1,PH2O)
    return T_SOL, PH2O_SOL

def sat_H2O_atmos(T):
    T0 = 273.15
    Pzero = 6.103E-3
    AMV = 18.
    VAPL = 597.3
    SUBL = 677.
    R = 1.9872
    A0 = 0.553
    BK = 1.38E-16
    
    HL = SUBL
    A = 0
    if T < T0:
        pass
    else:
        HL = VAPL + A0*T0
        A = A0

    P1 = Pzero * (T0/T)**(AMV*A/R)
    P2 = np.exp(AMV*HL/R * (1./T0 - 1./T))
    PV = P1 * P2 # bars
    return PV

def main():
    
    mu = 18.01534
    T_triple = 273.15
    T_critical = 647.0
    T_ref = 373.15 # K
    P_ref = 1.0142e6 # dynes/cm^2
    h = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref)

    T_data, P_data = get_data()
    guess = [3e10, 0, 2.4e10, 0]
    sol = h.optimize(guess, T_data, P_data)

    # save
    h.to_file('results/H2O_sat.yaml')

    entry = {}
    entry['T'] = np.linspace(200,600,100)
    entry['P'] = np.array([10.0**(6.8543 - 1897/(T-64.848)) for T in entry['T']])
    entry['label'] = 'Rimmer (2021)'
    data = [entry]

    entry = {}
    entry['T'] = np.linspace(T_data[0],642,100)
    entry['P'] = np.array([sat_H2O_atmos(T) for T in entry['T']])
    entry['label'] = 'Atmos code'
    data.append(entry)

    gas_g = 'H2O'
    fitting.plot_sat(h, gas_g, T_data, P_data, data)
    fitting.fit_thermo(gas_g, 100, 1000)

if __name__ == '__main__':
    main()