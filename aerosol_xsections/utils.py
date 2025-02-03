import h5py
import numpy as np
import miepython
import matplotlib.pyplot as plt

def write_file(filename, notes, wavelengths, radii, w0, qext, g0):

    with h5py.File('results/'+filename,'w') as f:

        dset = f.create_dataset("notes",(),dtype="S1000")
        dset[()] = '{:1000}'.format(notes).encode()

        dset = f.create_dataset("wavelengths", wavelengths.shape, 'f')
        dset[:] = wavelengths

        dset = f.create_dataset("radii", radii.shape, 'f')
        dset[:] = radii

        dset = f.create_dataset("w0", w0.T.shape, 'f')
        dset[:] = w0.T

        dset = f.create_dataset("qext", qext.T.shape, 'f')
        dset[:] = qext.T

        dset = f.create_dataset("g0", g0.T.shape, 'f')
        dset[:] = g0.T

def compute_mie_and_save(filename, notes, wavelength, m_real, m_imag, r_min, r_max, nrad, delete_fringe=True):

    radii, rup = get_r_grid_w_max(r_min, r_max, nrad)
    wavelength_nm = wavelength*1e3

    assert nrad == 50
    nw = len(wavelength)

    w0_all = np.zeros((nrad,nw),np.float64)
    qext_all = np.zeros((nrad,nw),np.float64)
    g_all = np.zeros((nrad,nw),np.float64)

    for i,radius in enumerate(radii):
        if i == 0:
            dr5 = (( rup[0] - radii[0] ) / 5.)
            rr = radii[0]
        else:
            dr5 = ( rup[i] - rup[i-1] ) / 5.
            rr = rup[i-1]

        if not delete_fringe:
            rr = radius
   
        for j in range(6):
            x = 2 * np.pi * rr / wavelength
            m = m_real - 1j*m_imag
            qext, qsca, qback, g = miepython.mie(m, x)
            w0 = qsca/qext
            
            w0_all[i,:] += w0/6
            qext_all[i,:] += qext/6
            g_all[i,:] += g/6
            rr += dr5
            if not delete_fringe:
                w0_all[i,:] *= 6
                qext_all[i,:] *= 6
                g_all[i,:] *= 6
                break

    write_file(filename, notes, wavelength_nm, radii, w0_all, qext_all, g_all)

def get_r_grid_w_max(r_min, r_max, n_radii):
    """
    Get spacing of radii to run Mie code

    r_min : float 
            Minimum radius to compute (cm)
    r_max : float 
            Maximum radius to compute (cm)
    n_radii : int
            Number of radii to compute 
    """

    radius = np.logspace(np.log10(r_min),np.log10(r_max),n_radii)
    rat = radius[1]/radius[0]
    rup = 2*rat / (rat+1) * radius

    return radius, rup

def save_plot(filename, radius):

    with h5py.File('results/'+filename,'r') as f:
        qext = f['qext'][:].T
        w0 = f['w0'][:].T
        g0 = f['g0'][:].T
        wavelength = f['wavelengths'][:]
        radii = f['radii'][:]

    ind = np.argmin(np.abs(radii-radius))
    fig, axs = plt.subplots(1,3,figsize=[11,3], sharex=False)
    axs[0].plot(wavelength, qext[ind,:], lw=0.5, c='k', marker='o', ms=1.5)
    axs[1].plot(wavelength, w0[ind,:], lw=0.5, c='k', marker='o', ms=1.5)
    axs[2].plot(wavelength, g0[ind,:], lw=0.5, c='k', marker='o', ms=1.5)

    for i,label in enumerate(['qext','w0','g0']):
        axs[i].set_ylabel(label)
    
    for ax in axs:
        ax.set_xlabel('Wavelength (nm)')
        ax.set_xscale('log')
        ax.grid(alpha=0.4)

    axs[0].text(0.0, 1.02, 'radius = %.4f microns'%(radii[ind]), \
            size = 12,ha='left', va='bottom',transform=axs[0].transAxes)
    
    plt.subplots_adjust(wspace=0.3)
    plt.savefig('figures/'+filename.strip('.h5')+'_%.4f.pdf'%(radii[ind]),bbox_inches='tight')

