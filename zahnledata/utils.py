
import numpy as np
from scipy import constants as const
from photochem.utils._format import yaml, flowmap, blockseqtrue
import decimal

class CustomDumper(yaml.Dumper):
    # modified from https://github.com/yaml/pyyaml/blob/master/lib3/yaml/representer.py#L171
    def represent_float(self, data):
        if data != data or (data == 0.0 and data == 1.0):
            value = '.nan'
        elif data == self.inf_value:
            value = '.inf'
        elif data == -self.inf_value:
            value = '-.inf'
        else:
            ctx = decimal.Context()
            ctx.prec = 7
            d1 = ctx.create_decimal(repr(float(data)))
            if np.abs(data) == 0:
                value = format(d1, 'f')
            elif np.abs(data) > 1e4:
                value = format(d1, 'e')
            elif np.abs(data) < 1e-3:
                value = format(d1, 'e')
            else:
                value = format(d1, 'f')
        return self.represent_scalar('tag:yaml.org,2002:float', value)

CustomDumper.add_representer(float, CustomDumper.represent_float)

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

def heat_capacity_shomate(coeffs, T):
    TT = T/1000.0
    cp = coeffs[0] + coeffs[1]*TT + coeffs[2]*TT**2 + \
         coeffs[3]*TT**3 + coeffs[4]/TT**2
    return cp # J/(mol K)

def heat_capacity_shomate_derivative(coeffs, T):
    TT = T/1000.0
    cp_dT = coeffs[1] + 2*coeffs[2]*TT + \
         3*coeffs[3]*TT**2 - 2*coeffs[4]/TT**3
    return cp_dT/1000

def gibbs_energy_shomate(coeffs, T):
    enthalpy = enthalpy_shomate(coeffs, T)
    entropy = entropy_shomate(coeffs, T)
    gibbs = enthalpy - T*entropy
    return gibbs

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

def heat_capacity_nasa9(coeffs, T):
    cp = (coeffs[0]*T**-2 + coeffs[1]*T**-1 + coeffs[2] + coeffs[3]*T \
         + coeffs[4]*T**2 + coeffs[5]*T**3 + coeffs[6]*T**4)*const.R
    return cp # J/mol/K

def gibbs_energy_nasa9(coeffs, T):
    enthalpy = enthalpy_nasa9(coeffs, T)
    entropy = entropy_nasa9(coeffs, T)
    gibbs = enthalpy - T*entropy
    return gibbs

def enthalpy_eval(thermo, T):
    ind = 0
    if T <= thermo['temperature-ranges'][0]:
        ind = 0
    
    elif T > thermo['temperature-ranges'][-1]:
        ind = -1

    else:
        for i in range(len(thermo['temperature-ranges'])-1):
            if T > thermo['temperature-ranges'][i] and T <= thermo['temperature-ranges'][i+1]:
                ind = i
                break
    
    if thermo['model'] == 'Shomate':
        enthalpy = enthalpy_shomate(thermo['data'][ind], T)
    elif thermo['model'] == 'NASA9':
        enthalpy = enthalpy_nasa9(thermo['data'][ind], T)
    else:
        raise Exception('Not found')

    return enthalpy # J/mol

def entropy_eval(thermo, T):

    ind = 0
    if T <= thermo['temperature-ranges'][0]:
        ind = 0
    
    elif T > thermo['temperature-ranges'][-1]:
        ind = -1

    else:
        for i in range(len(thermo['temperature-ranges'])-1):
            if T > thermo['temperature-ranges'][i] and T <= thermo['temperature-ranges'][i+1]:
                ind = i
                break
    
    if thermo['model'] == 'Shomate':
        entropy = entropy_shomate(thermo['data'][ind], T)
    elif thermo['model'] == 'NASA9':
        entropy = entropy_nasa9(thermo['data'][ind], T)
    else:
        raise Exception('Not found')

    return entropy # J/mol/K

def heat_capacity_eval(thermo, T):
    ind = 0
    if T <= thermo['temperature-ranges'][0]:
        ind = 0
    
    elif T > thermo['temperature-ranges'][-1]:
        ind = -1

    else:
        for i in range(len(thermo['temperature-ranges'])-1):
            if T > thermo['temperature-ranges'][i] and T <= thermo['temperature-ranges'][i+1]:
                ind = i
                break

    if thermo['model'] == 'Shomate':
        heat_capacity = heat_capacity_shomate(thermo['data'][ind], T)
    elif thermo['model'] == 'NASA9':
        heat_capacity = heat_capacity_nasa9(thermo['data'][ind], T)
    else:
        raise Exception('Not found')

    return heat_capacity # J/mol/K
    
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
        
        self.mech = {}
        self.mech['species'] = species
                
    def enthalpy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return enthalpy_eval(thermo, T)
    
    def entropy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return entropy_eval(thermo, T)
    
    def gibbs_energy(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return gibbs_energy_eval(thermo, T)
    
    def heat_capacity(self, sp, T):
        thermo = self.mech['species'][sp]['thermo']
        return heat_capacity_eval(thermo, T)
