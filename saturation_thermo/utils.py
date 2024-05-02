
import numpy as np
from scipy import constants as const
from scipy import optimize
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
            d1 = ctx.create_decimal(repr(data))
            if np.abs(data) > 1e4:
                value = format(d1, 'e')
            else:
                value = format(d1, 'f')
        return self.represent_scalar('tag:yaml.org,2002:float', value)

CustomDumper.add_representer(float, CustomDumper.represent_float)

class SaturationPropertiesFitter():
    R = 8.31446261815324e7 # ideal gas in cgs units
    
    def __init__(self, mu, T_triple, T_critical, T_ref, P_ref):
        self.mu = mu # molar weight
        self.T_triple = T_triple
        self.T_critical = T_critical
        self.T_ref = T_ref
        if T_ref < T_triple:
            raise Exception("T_ref must be bigger than T_triple")
        self.P_ref = P_ref
    
    def _latent_heat(self, A, B, T):
        return A + B*T
    
    def latent_heat(self, A_v, B_v, A_s, B_s, A_c, B_c, T):
        if T > self.T_critical:
            return self._latent_heat(A_c, B_c, T)
        elif self.T_critical >= T > self.T_triple:
            return self._latent_heat(A_v, B_v, T)
        elif T <= self.T_triple:
            return self._latent_heat(A_s, B_s, T)
    
    def _integral(self, A, B, T):
        return -A/T + B*np.log(T)
    
    def sat_pressure(self, A_v, B_v, A_s, B_s, A_c, B_c, T):
        
        if T > self.T_critical:
            tmp = (self._integral(A_v, B_v, self.T_critical) - self._integral(A_v, B_v, self.T_ref)) + \
                  (self._integral(A_c, B_c, T) - self._integral(A_c, B_c, self.T_critical))
        elif self.T_critical >= T > self.T_triple:
            tmp = self._integral(A_v, B_v, T) - self._integral(A_v, B_v, self.T_ref)
        elif T <= self.T_triple:
            tmp = (self._integral(A_v, B_v, self.T_triple) - self._integral(A_v, B_v, self.T_ref)) + \
                  (self._integral(A_s, B_s, T) - self._integral(A_s, B_s, self.T_triple))
            
        Ps = self.P_ref*np.exp((self.mu/self.R)*(tmp))
        return Ps
    
    def objective(self, x, T_sol, P_sol):
        A_v,B_v,A_s,B_s = x   
        tmp = np.array([self.sat_pressure(A_v,B_v,A_s,B_s,0.0,0.0,T) for T in T_sol])
        out = np.log10(tmp) - np.log10(P_sol)
        return out
    
    def optimize(self, guess, T_sol, P_sol):
        sol = optimize.root(self.objective, guess, method='lm', args=(T_sol,P_sol,), options={'maxiter':100000})
        if not sol.success:
            raise Exception('root solve failed')
            
        self.A_v = sol.x[0]
        self.B_v = sol.x[1]
        self.A_s = sol.x[2]
        self.B_s = sol.x[3]

        # Deal with supercritical
        P_c = self.sat_pressure(self.A_v,self.B_v,self.A_s,self.B_s,0,0,self.T_critical)
        T = self.T_critical + 5
        P = 10.0**(np.log10(P_c) + 2)
        self.A_c = (self.R/self.mu)*(1/self.T_critical - 1/T)**(-1)*np.log(P/P_c)
        self.B_c = 0.0

        return sol
    
    def to_dict(self):
        sol = {}
        sol['model'] = 'LinearLatentHeat'
        sol['mu'] = float(self.mu)
        sol['T-ref'] = float(self.T_ref)
        sol['P-ref'] = float(self.P_ref)
        sol['T-triple'] = float(self.T_triple)
        sol['T-critical'] = float(self.T_critical)
        sol['vaporization'] = {}
        sol['vaporization']['A'] = float(self.A_v)
        sol['vaporization']['B'] = float(self.B_v)
        sol['sublimation'] = {}
        sol['sublimation']['A'] = float(self.A_s)
        sol['sublimation']['B'] = float(self.B_s)
        sol['super-critical'] = {}
        sol['super-critical']['A'] = float(self.A_c)
        sol['super-critical']['B'] = float(self.B_c)
        return sol
    
    def to_file(self, filename):
        sol = self.to_dict()
        sol['vaporization'] = flowmap(sol['vaporization'])
        sol['sublimation'] = flowmap(sol['sublimation'])
        sol['super-critical'] = flowmap(sol['super-critical'])
        with open(filename,'w') as f:
            yaml.dump(sol,f,Dumper=CustomDumper,sort_keys=False,width=70)

class SaturationProperties():

    R = 8.31446261815324e7 # ideal gas in cgs units

    def __init__(self):
        pass
    
    def init1(self, mu, T_triple, T_critical, T_ref, P_ref, consts):
        self.mu = mu # molar weight
        self.T_triple = T_triple
        self.T_critical = T_critical
        self.T_ref = T_ref
        if T_ref < T_triple:
            raise Exception("T_ref must be bigger than T_triple")
        self.P_ref = P_ref
        self.A_v = consts[0]
        self.B_v = consts[1]
        self.A_s = consts[2]
        self.B_s = consts[3]
        self.A_c = consts[4]
        self.B_c = consts[5]

    def init2(self, filename):
        with open(filename,'r') as f:
            tmp = yaml.load(f,Loader=yaml.Loader)
        self.mu = tmp['mu']
        self.T_triple = tmp['T-triple']
        self.T_critical = tmp['T-critical']
        self.T_ref = tmp['T-ref']
        if self.T_ref < self.T_triple:
            raise Exception("T_ref must be bigger than T_triple")
        self.P_ref = tmp['P-ref']
        self.A_v = tmp['vaporization']['A']
        self.B_v = tmp['vaporization']['B']
        self.A_s = tmp['sublimation']['A']
        self.B_s = tmp['sublimation']['B']
        self.A_c = tmp['super-critical']['A']
        self.B_c = tmp['super-critical']['B']
    
    def _latent_heat(self, A, B, T):
        return A + B*T
    
    def latent_heat(self, T):
        if T > self.T_critical:
            return self._latent_heat(self.A_c, self.B_c, T)
        elif self.T_critical >= T > self.T_triple:
            return self._latent_heat(self.A_v, self.B_v, T)
        elif T <= self.T_triple:
            return self._latent_heat(self.A_s, self.B_s, T)
    
    def _integral(self, A, B, T):
        return -A/T + B*np.log(T)
    
    def sat_pressure(self, T):
        
        if T > self.T_critical:
            tmp = (self._integral(self.A_v, self.B_v, self.T_critical) - self._integral(self.A_v, self.B_v, self.T_ref)) + \
                  (self._integral(self.A_c, self.B_c, T) - self._integral(self.A_c, self.B_c, self.T_critical))
        elif self.T_critical >= T > self.T_triple:
            tmp = self._integral(self.A_v, self.B_v, T) - self._integral(self.A_v, self.B_v, self.T_ref)
        elif T <= self.T_triple:
            tmp = (self._integral(self.A_v, self.B_v, self.T_triple) - self._integral(self.A_v, self.B_v, self.T_ref)) + \
                  (self._integral(self.A_s, self.B_s, T) - self._integral(self.A_s, self.B_s, self.T_triple))
            
        Ps = self.P_ref*np.exp((self.mu/self.R)*(tmp))
        return Ps

def create_CO2(): 
    mu_CO2 = 44.01
    T_triple = 216.58
    T_critical = 304.13
    T_ref = 250.0 # K
    P_ref = 17843676.678142548 # dynes/cm^2
    consts = [4.656475e+09, -3.393595e+06, 6.564668e+09, -3.892217e+06, 1.635908e+11, 0.0]
    sat_CO2 = SaturationProperties()
    sat_CO2.init1(mu_CO2, T_triple, T_critical, T_ref, P_ref, consts)
    return sat_CO2

def create_H2O(): 
    mu_H2O = 18.01534
    T_triple = 273.15
    T_critical = 647.0
    T_ref = 373.15 # K
    P_ref = 1.0142e6 # dynes/cm^2
    consts = [2.841421e+10, -1.399732e+07, 2.746884e+10, 4.181527e+06, 1.793161e+12, 0.0]
    sat_H2O = SaturationProperties()
    sat_H2O.init1(mu_H2O, T_triple, T_critical, T_ref, P_ref, consts)
    return sat_H2O

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

def Psat_eqn(Gs, Gg, T):
    DGr = Gs - Gg
    Keq = np.exp(-DGr/(const.R*T))
    Psat = 1/Keq
    return Psat # bar

def latent_heat_eqn(Hs, Hg):
    return Hg - Hs

def Psat_eval(thermo_s, thermo_g, T):
    Gs = gibbs_energy_eval(thermo_s, T)
    Gg = gibbs_energy_eval(thermo_g, T)
    return Psat_eqn(Gs, Gg, T)

def latent_heat_eval(thermo_s, thermo_g, T):
    Hs = enthalpy_eval(thermo_s, T)
    Hg = enthalpy_eval(thermo_g, T)
    return latent_heat_eqn(Hs, Hg)

def saturation_lodders(a, b, T):
    return 10.0**(a + b/T)*1e6 # dynes/cm^2

def saturation_antoine(A, B, C, T):
    return 10.0**(A - (B/(T+C)))*1e6 # dynes/cm^2

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
    
    def Psat(self, sp_s, sp_g, T):
        Gs = self.gibbs_energy(sp_s,T)
        Gg = self.gibbs_energy(sp_g,T)
        return Psat_eqn(Gs, Gg, T)
    
    def latent_heat(self, sp_s, sp_g, T):
        Hs = self.enthalpy(sp_s,T)
        Hg = self.enthalpy(sp_g,T)
        return latent_heat_eqn(Hs, Hg)

    
def Format_species(sp):
            
    order = ['name', 'composition', 'condensate', 'triple-temperature', 'critical-temperature', 'thermo','note']
    copy = sp.copy()
    sp.clear()
    for key in order:
        if key in copy.keys():
            sp[key] = copy[key]
                
    sp['composition'] = flowmap(sp['composition'])
    if 'thermo' in sp.keys():
        
        order = ['model', 'reference-pressure','temperature-ranges','data']
        copy = sp['thermo'].copy()
        sp['thermo'].clear()
        for key in order:
            if key in copy.keys():
                sp['thermo'][key] = copy[key]
            
        sp['thermo']['temperature-ranges'] = blockseqtrue(sp['thermo']['temperature-ranges'])
        
        sp['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(sp['thermo']['data'])]

    return sp