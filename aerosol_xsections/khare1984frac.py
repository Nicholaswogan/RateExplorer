import numpy as np
import utils
import fractal

def main():
    _, wavelength, m_real, m_imag = np.loadtxt('data/khare_tholins.dat',skiprows=13).T
    wavelength_save = np.logspace(np.log10(0.0207), np.log10(920.0),500) # microns
    m_real_save = np.interp(wavelength_save, wavelength, m_real)
    m_imag_save = np.interp(wavelength_save, wavelength, m_imag)

    r_min = 0.05
    r_max = 5
    nrad = 50
    a = 1
    df = 2
    rmon = 0.05
    nprocess = 4
    filename = 'mie_khare1984frac.h5'
    notes = """Mie optical properties for hydrocarbon aerosols based on Khare et al. (1984).

wavelengths: Wavelength [nm]
radii: Particle radii [um]
w0: Single scattering albedo [unitless], dimensions (len(radii),len(wavelengths))
qext: Extinction [1/particle], dimensions (len(radii),len(wavelengths))
g0: Asymmetry factor [unitless], dimensions (len(radii),len(wavelengths))
"""
    fractal.compute_frac_and_save(filename, notes, wavelength_save, m_real_save, m_imag_save, 
                                  r_min, r_max, nrad, 
                                  a, df, rmon, nprocess)

    utils.save_plot(filename, 0.1)
    utils.save_plot(filename, 0.5)
    utils.save_plot(filename, 1)
    utils.save_plot(filename, 5)
    

if __name__ == '__main__':
    main()