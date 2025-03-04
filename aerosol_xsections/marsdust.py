import numpy as np
import utils

def main():
    wavelength, m_real, m_imag = np.loadtxt('data/mars043i.dat',skiprows=3).T
    wavelength_save = np.logspace(np.log10(np.min(wavelength)), np.log10(np.max(wavelength)),500) # microns
    m_real_save = np.interp(wavelength_save, wavelength, m_real)
    m_imag_save = np.interp(wavelength_save, wavelength, m_imag)

    r_min = 0.05
    r_max = 50
    nrad = 50
    filename = 'marsdust.h5'
    notes = """Mie optical properties for Mars dust based on Wolff et al. (2009).

wavelengths: Wavelength [nm]
radii: Particle radii [um]
w0: Single scattering albedo [unitless], dimensions (len(radii),len(wavelengths))
qext: Extinction [1/particle], dimensions (len(radii),len(wavelengths))
g0: Asymmetry factor [unitless], dimensions (len(radii),len(wavelengths))
"""
    utils.compute_mie_and_save(filename, notes, wavelength_save, m_real_save, m_imag_save, r_min, r_max, nrad, delete_fringe=False)

    utils.save_plot(filename, 0.1)
    utils.save_plot(filename, 0.5)
    utils.save_plot(filename, 1)
    utils.save_plot(filename, 10)
    

if __name__ == '__main__':
    main()