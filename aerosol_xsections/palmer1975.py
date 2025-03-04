import numpy as np
import utils

def main():

    # Read
    with open('data/palmer_williams_h2so4.dat','r') as fil:
        lines = fil.readlines()
    wavelength = [] # microns
    m_real = []
    m_imag = []
    for line in lines[16:243]:
        tmp = line.split()
        wavelength.append(float(tmp[1]))
        m_real.append(float(tmp[7]))
    for line in lines[248:-1]:
        tmp = line.split()
        m_imag.append(float(tmp[7]))
    wavelength = np.array(wavelength) # this is micro meter
    m_real = np.array(m_real)
    m_imag = np.array(m_imag)

    # Sort
    inds = np.argsort(wavelength)
    wavelength = wavelength[inds].copy()
    m_real = m_real[inds].copy()
    m_imag = m_imag[inds].copy()

    # Interpolate to new grid
    wavelength_save = np.logspace(np.log10(1e-3), np.log10(1000.0), 500) # microns
    m_real_save = np.interp(wavelength_save, wavelength, m_real)
    m_imag_save = np.interp(wavelength_save, wavelength, m_imag)

    r_min = 0.01
    r_max = 100.0
    nrad = 50
    filename = 'palmer1975.h5'
    notes = """Mie optical properties for sulfuric acid aerosols based on Palmer and Williams (1975).

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
    

if __name__ == '__main__':
    main()