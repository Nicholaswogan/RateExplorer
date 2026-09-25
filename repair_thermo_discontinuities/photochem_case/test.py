import planets

from photochem.utils import stars
from photochem.extensions import gasgiants as gasgiants_0_9_0
from photochem._clima import rebin, rebin_with_errors
import gasgiants_0_6_7 as gasgiants_0_6_7

import numpy as np
from scipy import interpolate
from astropy import constants
from scipy import constants as const

import astropy.units as u
import pickle
import os
import re
from copy import deepcopy

def get_comp_yang(metallicity, x_acc=None, CtoO=None, N_depletion=None):
    """Yang and Hu (2024) method for computing atomic composition.
    x_acc is [H2]/([H2] + [H2O]). So 
    """

    if N_depletion is not None and x_acc is not None:
        raise Exception('')
    if CtoO is not None and x_acc is not None:
        raise Exception('')
    
    comp = {
        'O': 0.000495,
        'C': 0.000272,
        'N': 0.000065,
        'S': 0.000013,
        'He': 0.077379,
        'H': 0.921775
    }
    tot = sum(comp.values())
    for key in comp:
        comp[key] /= tot
    for key in comp:
        if key not in ['H','He']:
            comp[key] *= metallicity
    tot = sum(comp.values())
    for key in comp:
        comp[key] /= tot

    # Acreete H2O
    if x_acc is not None:
        a = comp['H'] + comp['He'] + comp['O']
        comp['H'] = (2*a)/(2 + (1 - x_acc) + 2/11.91)
        comp['He'] = (2*a)/((2 + (1 - x_acc) + 2/11.91)*11.91)
        comp['O'] = a*(1 - x_acc)/(2 + (1 - x_acc) + 2/11.91)

    # Apply C/O ratio
    if CtoO is not None:
        x = CtoO*(comp['C']/comp['O'])
        a = (x*comp['O'] - comp['C'])/(1 + x)
        comp['C'] = comp['C'] + a
        comp['O'] = comp['O'] - a

    # Deplete N
    if N_depletion is not None:
        comp['N'] /= N_depletion
        
    tot = sum(comp.values())
    for key in comp:
        comp[key] /= tot

    return comp

def initialize_photochem(spectrum, planet_mass, planet_radius, metallicity, x_acc, CtoO, N_depletion, climate_filename, Kzz, initial_cond_with_quenching):

    gasgiants = gasgiants_0_9_0

    pc = gasgiants.EvoAtmosphereGasGiant(
        'photochem_rxns.yaml',
        spectrum,
        planet_mass*constants.M_earth.value*1e3, # grams
        planet_radius*constants.R_earth.value*1e2, # cm
        solar_zenith_angle=60,
        thermo_file='photochem_thermo.yaml'
    )
    pc.gdat.verbose = True
    pc.var.verbose = 1
    # pc.var.jacobian_method = 1
    pc.gdat.TOA_pressure_avg = 0.01
    pc.var.nsteps_before_giveup = 8000
    pc.gdat.max_total_step = 8000
    # print(pc.gdat.max_total_step)

    # Set particle radius
    particle_radius = pc.var.particle_radius
    particle_radius[:,:] = 1e-3
    pc.var.particle_radius = particle_radius
    # pc.update_vertical_grid(TOA_alt=pc.var.top_atmos)

    # Condensation and Evaporation parameters
    for i in range(len(pc.var.cond_params)):
        pc.var.cond_params[i].smooth_factor = 2
        pc.var.cond_params[i].k_cond = 1000
        pc.var.cond_params[i].k_evap = 0

    # Set composition
    comp = get_comp_yang(metallicity, x_acc=x_acc, CtoO=CtoO, N_depletion=N_depletion)
    molfracs_atoms = np.empty(len(pc.gdat.gas.atoms_names))
    for i,atom in enumerate(pc.gdat.gas.atoms_names):
        molfracs_atoms[i] = comp[atom]
    pc.gdat.gas.molfracs_atoms_sun = molfracs_atoms

    with open(climate_filename,'rb') as f:
        out = pickle.load(f)
    P = out['pressure'][::-1].copy()*1e6
    T = out['temperature'][::-1].copy()
    if np.max(P) <= 1e10:
        P_append = np.arange(np.log10(np.max(P)), 10.01,.1)[::-1][:-1]
        T1 = interpolate.interp1d(np.log10(P[::-1]), T[::-1], fill_value='extrapolate')(P_append)
        P = np.append(10.0**P_append,P)
        T = np.append(T1,T)
    Kzz1 = np.ones(P.shape[0])*Kzz

    pc.gdat.initial_cond_with_quenching = initial_cond_with_quenching

    pc.initialize_to_climate_equilibrium_PT(P, T, Kzz1, 1, 1)

    return pc

def run():

    pl = planets.TOI1231b
    name = 'TOI1231b'
    pc = initialize_photochem(
        spectrum=f'{name}_spectrum.txt',
        planet_mass=pl.mass,
        planet_radius=pl.radius,
        metallicity=80,
        x_acc=None,
        CtoO=1.0,
        N_depletion=None,
        climate_filename=f'{name}_MH=2.000_CO=1.000_Tint=50.0.pkl',
        Kzz=1e7,
        initial_cond_with_quenching=True
    )

    converged = pc.find_steady_state()
    print(converged)

if __name__ == '__main__':
    run()