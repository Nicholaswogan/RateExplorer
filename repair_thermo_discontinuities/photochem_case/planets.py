import warnings
warnings.filterwarnings("ignore")
import numpy as np
from photochem.utils import stars

class Star:
    radius : float # relative to the sun
    Teff : float # K
    metal : float # log10(M/H)
    kmag : float
    logg : float
    planets : dict # dictionary of planet objects

    def __init__(self, radius, Teff, metal, kmag, logg, planets):
        self.radius = radius
        self.Teff = Teff
        self.metal = metal
        self.kmag = kmag
        self.logg = logg
        self.planets = planets
        
class Planet:
    radius : float # in Earth radii
    mass : float # in Earth masses
    Teq : float # Equilibrium T in K
    transit_duration : float # in seconds
    a: float # semi-major axis in AU
    a: float # semi-major axis in AU
    stellar_flux: float # W/m^2
    
    def __init__(self, radius, mass, Teq, transit_duration, a, stellar_flux):
        self.radius = radius
        self.mass = mass
        self.Teq = Teq
        self.transit_duration = transit_duration
        self.a = a
        self.stellar_flux = stellar_flux

# Greklek-McKeon et al. 2024

TOI1266b = Planet(
    radius=2.52,
    mass=4.46,
    Teq=415,
    transit_duration=2.1288*60*60, # Exo.Mast
    a=0.0730,
    stellar_flux=stars.equilibrium_temperature_inverse(415, 0.0)
)

TOI1266c = Planet(
    radius=1.98,
    mass=3.17,
    Teq=346,
    transit_duration=1.951*60*60, # Exo.Mast
    a=0.1050,
    stellar_flux=stars.equilibrium_temperature_inverse(346, 0.0)
)

TOI1266 = Star(
    radius=0.4232,
    Teff=3563,
    metal=-0.109,
    kmag=8.84,
    logg=4.826,
    planets={'b': TOI1266b, 'c': TOI1266c}
)


K2_18b = Planet(
    radius=2.610, # Benneke et al. (2019)
    mass=8.63, # Cloutier et al. (2019)
    Teq=278.7, # Matches stellar constant
    transit_duration=2.682*60*60, # Exo.Mast
    a=0.15910, # Benneke et al. (2019)
    stellar_flux=1368.0 # Benneke et al. (2019)
)

K2_18 = Star(
    radius=0.4445, # Benneke et al. (2019)
    Teff=3457, # Benneke et al. (2017)
    metal=0.12, # Exo.Mast
    kmag=8.9, # Exo.Mast
    logg=4.79, # Exo.Mast
    planets={'b':K2_18b}
)

# Piaulet-Ghorayeb et al. 2024

GJ9827d = Planet(
    radius=1.98,
    mass=3.02,
    Teq=675,
    transit_duration=1.223, # Rice et al. (2019)
    a=0.053,
    stellar_flux=stars.equilibrium_temperature_inverse(675, 0.0)
)

GJ9827 = Star(
    radius=0.58,
    Teff=4236,
    metal=-0.29,
    kmag=7.2, # exo.Mast
    logg=4.70,
    planets={'d': GJ9827d}
)

# Van Eylen et al. 2021

TOI270d = Planet(
    radius=2.133,
    mass=4.78,
    Teq=387,
    transit_duration=2.148*60*60, # Gunther et al. (2019)
    a=0.07210,
    stellar_flux=stars.equilibrium_temperature_inverse(387, 0.0)
)

TOI270 = Star(
    radius=0.378,
    Teff=3506,
    metal=-0.20,
    kmag=8.3, # exo.Mast
    logg=4.872,
    planets={'b': TOI270d}
)


# NASA Exoplanet Archive PSCompPars (accessed 2026-09-24).
# The discovery paper linked by NEA reports a 3-sigma upper mass limit of
# 9.9 M_earth for TOI-6478 b. NEA's 19.2 M_earth composite value is instead
# a calculated mass-radius estimate, so the measured upper limit is used here.
TOI6478b = Planet(
    radius=4.60,
    mass=9.9,  # 3-sigma upper limit; Scott et al. (2025)
    Teq=223.5,  # Zero-albedo value consistent with the stellar properties and a
    transit_duration=2.900 * 60 * 60,
    a=0.1136,
    stellar_flux=stars.equilibrium_temperature_inverse(223.5, 0.0)
)

TOI6478 = Star(
    radius=0.234,
    Teff=3230,
    metal=-0.53,
    kmag=10.657,
    logg=5.016,
    planets={'b': TOI6478b}
)

# NASA Exoplanet Archive PSCompPars (accessed 2026-09-24).
TOI1231b = Planet(
    radius=3.65,
    mass=15.4,
    Teq=329.6,
    transit_duration=3.240 * 60 * 60,
    a=0.1288,
    stellar_flux=stars.equilibrium_temperature_inverse(329.6, 0.0)
)

TOI1231 = Star(
    radius=0.476,
    Teff=3553,
    metal=0.041,
    kmag=8.069,
    logg=4.767,
    planets={'b': TOI1231b}
)
