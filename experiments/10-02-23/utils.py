
import yaml
import numpy as np
from scipy import constants as const

def arrhenius_rate(A, b, Ea, T):
    return A*T**b*np.exp(-Ea/T)

def falloff_rate(kinf, Pr):
    k = kinf * (Pr / (1.0 + Pr))
    return k

def enthalpy_shomate(coeffs, T):
    TT = T/1000.0
    enthalpy = coeffs[0]*TT + (coeffs[1]*TT**2)/2.0 \
             + (coeffs[2]*TT**3)/3.0  + (coeffs[3]*TT**4)/4.0 \
             - coeffs[4]/TT + coeffs[5]
    return enthalpy*1000.0 # J/mol

def entropy_shomate(coeffs, T):
    TT = T/1000.0
    entropy = coeffs[0]*np.log(TT) + coeffs[1]*TT \
            + (coeffs[2]*TT**2)/2.0 + (coeffs[3]*TT**3)/3.0 \
            - coeffs[4]/(2.0 * TT**2) + coeffs[6]
    return entropy # J/mol/K

def enthalpy_nasa9(coeffs, T):
    enthalpy = (- coeffs[0]*T**(-2.0) + coeffs[1]*np.log(T)/T \
                + coeffs[2] + coeffs[3]*T/2.0 + coeffs[4]*T**(2.0)/3.0 \
                + coeffs[5]*T**(3.0)/4.0 + coeffs[6]*T**(4.0)/5.0 \
                + coeffs[7]/T)*T*const.R
    return enthalpy # J/mol

def entropy_nasa9(coeffs, T):
    entropy = (- coeffs[0]*T**(-2.0)/2.0 - coeffs[1]*T**(-1.0) \
               + coeffs[2]*np.log(T) + coeffs[3]*T + coeffs[4]*T**(2.0)/2.0 \
               + coeffs[5]*T**(3.0)/3.0 + coeffs[6]*T**(4.0)/4.0 \
               + coeffs[8])*const.R
    return entropy # J/mol/K

def enthalpy_eval(thermo, T):
    found = False
    for i in range(len(thermo['temperature-ranges'])-1):
        if T > thermo['temperature-ranges'][i] and T <= thermo['temperature-ranges'][i+1]:
            found = True
            if thermo['model'] == 'Shomate':
                enthalpy = enthalpy_shomate(thermo['data'][i], T)
            elif thermo['model'] == 'NASA9':
                enthalpy = enthalpy_nasa9(thermo['data'][i], T)
                pass
            else:
                raise Exception('Not found')

    if not found:
        raise Exception('Not found')

    return enthalpy # J/mol

def entropy_eval(thermo, T):
    found = False
    for i in range(len(thermo['temperature-ranges'])-1):
        if T > thermo['temperature-ranges'][i] and T <= thermo['temperature-ranges'][i+1]:
            found = True
            if thermo['model'] == 'Shomate':
                entropy = entropy_shomate(thermo['data'][i], T)
            elif thermo['model'] == 'NASA9':
                entropy = entropy_nasa9(thermo['data'][i], T)
                pass
            else:
                raise Exception('Not found')
            
    if not found:
        raise Exception('Not found')

    return entropy # J/mol/K
    
def gibbs_energy_eval(thermo, T):
    enthalpy = enthalpy_eval(thermo, T)
    entropy = entropy_eval(thermo, T)
    gibbs = enthalpy - T*entropy
    return gibbs # J/mol

class ReactionExplorer():
    
    def __init__(self, mech_file):
        with open(mech_file,'r') as f:
            data = yaml.load(f,Loader=yaml.Loader)
        
        species = {}
        for i in range(len(data['species'])):
            sp = data['species'][i]['name']
            species[sp] = data['species'][i]
        
        reactions = []
        for i in range(len(data['reactions'])):
            eq = data['reactions'][i]['equation']
            eq = eq.replace('(','').replace(')','')
            react, prod = eq.split('<=>')
            react = react.split('+')
            react = [a.strip() for a in react]
            if 'M' in react:
                react.remove('M')
            prod = prod.split('+')
            prod = [a.strip() for a in prod]
            if 'M' in prod:
                prod.remove('M')
            rx = data['reactions'][i]
            if 'type' not in rx:
                rx['type'] = 'elementary'

            rx['react'] = react
            rx['prod'] = prod
            reactions.append(rx)
        
        self.mech = {}
        self.mech['species'] = species
        self.mech['reactions'] = reactions

        if data['species'][0]['thermo']['model'] == 'Shomate':
            self.reference_pressure = 1.0e6
        elif data['species'][0]['thermo']['model'] == 'NASA9':
            self.reference_pressure = 1.01325e6

    def evaluate_rates(self, T, P):
        '''P in bars'''
        # compute the density
        k_boltz = const.k*1e7
        den = 1e6*P/(k_boltz*T)

        rates = []
        for i,rx in enumerate(self.mech['reactions']):
            if rx['type'] == 'elementary':
                rates.append(None)
            elif rx['type'] == 'three-body':
                rates.append(None)
            elif rx['type'] == 'falloff':
                tmp = rx['low-P-rate-constant']
                k0_f = arrhenius_rate(tmp['A'], tmp['b'], tmp['Ea'], T)
                tmp = rx['high-P-rate-constant']
                kinf_f = arrhenius_rate(tmp['A'], tmp['b'], tmp['Ea'], T)
                Pr = k0_f*den/kinf_f
                # k_f = falloff_rate(kinf_f, Pr)
                k_f = den*k0_f/(1+(k0_f*den/kinf_f))

                # compute gibbs energy
                gibbR_forward = 0.0
                for i,sp in enumerate(rx['react']):
                    thermo = self.mech['species'][sp]['thermo']
                    gibbR_forward += gibbs_energy_eval(thermo, T)

                gibbP_forward = 0.0
                for i,sp in enumerate(rx['prod']):
                    thermo = self.mech['species'][sp]['thermo']
                    gibbP_forward += gibbs_energy_eval(thermo, T)

                Dg_forward = gibbP_forward - gibbR_forward
                l = len(rx['react'])
                m = len(rx['prod'])
                reverse_factor = (1.0/np.exp(-Dg_forward/(const.R * T))) * (k_boltz*T/self.reference_pressure)**(m-l)

                out = {}
                out['equation'] = rx['equation']
                out['k0_f'] = k0_f
                out['kinf_f'] = kinf_f
                out['k_f'] = k_f
                out['k0_r'] = k0_f*reverse_factor
                out['kinf_r'] = kinf_f*reverse_factor
                out['k_r'] = k_f*reverse_factor
                rates.append(out)

        return rates
                
    def enthalpy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return enthalpy_eval(thermo, T)
    
    def entropy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return entropy_eval(thermo, T)
    
    def gibbs_energy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return gibbs_energy_eval(thermo, T)
    
