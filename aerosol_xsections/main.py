import numpy as np
import utils

def read_khare1984():
    _, wavelength, m_real, m_imag = np.loadtxt('data/khare_tholins.dat',skiprows=13).T
    return wavelength, m_real, m_imag

def read_marsdust():
    wavelength, m_real, m_imag = np.loadtxt('data/mars043i.dat',skiprows=3).T
    return wavelength, m_real, m_imag

def read_palmer1975():
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

    return wavelength, m_real, m_imag

def make_metadata(description, citation, wavelengths):
    return {
        'description': f'Mie optical properties for {description}',
        'citation': citation,
        'source_wavelength_min_um': float(wavelengths.min()),
        'source_wavelength_max_um': float(wavelengths.max()),
        'refractive_index_extrapolation': (
            'Hold real and imaginary refractive indices at their nearest '
            'source values outside the source wavelength range, then compute '
            'Mie optical properties.'
        ),
    }

def main():

    wavelength_save = np.logspace(np.log10(0.01), np.log10(1000.0), 1000)
    r_min = 0.01
    r_max = 500.0
    nrad = 80

    inputs = {
        'khare1984': {
            'description': 'hydrocarbon aerosols',
            'citation': 'Khare et al. (1984)',
            'data': read_khare1984(),
            'plot_radii': (0.1, 0.5, 1),
        },
        'marsdust': {
            'description': 'Mars dust',
            'citation': 'Wolff et al. (2009)',
            'data': read_marsdust(),
            'plot_radii': (0.1, 0.5, 1, 10),
        },
        'palmer1975': {
            'description': 'sulfuric acid aerosols',
            'citation': 'Palmer and Williams (1975)',
            'data': read_palmer1975(),
            'plot_radii': (0.1, 0.5, 1),
        },
    }

    for name, entry in inputs.items():
        wavelength, m_real, m_imag = entry['data']
        m_real_save = np.interp(wavelength_save, wavelength, m_real)
        m_imag_save = np.interp(wavelength_save, wavelength, m_imag)
        filename = f'mie_{name}.h5'
        metadata = make_metadata(entry['description'], entry['citation'], wavelength)

        utils.compute_mie_and_save(
            filename, metadata, wavelength_save, m_real_save, m_imag_save,
            r_min, r_max, nrad, delete_fringe=False,
        )
        for radius in entry['plot_radii']:
            utils.save_plot(filename, radius)

if __name__ == '__main__':
    main()
