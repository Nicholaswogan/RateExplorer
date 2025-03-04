import numpy as np
import utils
import ctypes as ct
import os
import multiprocessing as mp
from tqdm import tqdm

def number_of_monomers(R, R_mon):
    return (R/R_mon)**3

def make_fractal_meanfield():
    rootdir = os.path.dirname(os.path.realpath(__file__))+'/'
    files = os.listdir('./')
    filename = [a for a in files if 'libfractal' in a][0]
    libfractal = ct.CDLL(rootdir+filename)
    fractal_meanfield = libfractal.fractal_meanfield
    fractal_meanfield.argtypes = [
        ct.c_bool,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.c_double,
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int)
    ]
    fractal_meanfield.restype = None
    return fractal_meanfield
    
fractal_meanfield_wrapper = make_fractal_meanfield()

def fractal_meanfield_single(wv, k, n, a, df, r, rmon):

    nmon = number_of_monomers(r, rmon)

    qext = ct.c_double(0.0)
    qsca = ct.c_double(0.0)
    gfac = ct.c_double(0.0)
    rc = ct.c_int(0)
    
    fractal_meanfield_wrapper(
        False,
        wv,
        k,
        n,
        nmon,
        a,
        df,
        rmon,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ct.byref(qext),
        ct.byref(qsca),
        ct.byref(gfac),
        ct.byref(rc)
    )
    if rc.value < 0:
        return np.nan, np.nan, np.nan
    w0 = qsca.value/qext.value
    return qext.value, w0, gfac.value

def fractal_meanfield(wv, k, n, a, df, r, rmon, nprocess=None):
    """Fractal mean field.

    Parameters
    ----------
    wv : ndarray[size=1, double]
        Wavelength in microns.
    k : ndarray[size=1, double]
        Imaginary index of refraction
    n : ndarray[size=1, double]
        Real index of of refraction
    a : float
        Packing coefficient
    df : float
        Fractal dimension
    r : float
        Particle radius (um)
    rmon : float
        Monomer radius (um)

    Returns
    -------
    qext : ndarray[size=1, double]
        Extinction coefficient
    w0 : ndarray[size=1, double]
        Single scattering albedo
    g0 : ndarray[size=1, double]
        Asymmetry factor
    """
    if nprocess is None:
        nprocess = mp.cpu_count()
    with mp.Pool(nprocess) as pool:
        jobs = []
        for i in range(len(wv)):
            job = pool.apply_async(fractal_meanfield_single, (wv[i], k[i], n[i], a, df, r, rmon))
            jobs.append(job)
        res = []
        for job in tqdm(jobs):
            res.append(job.get())
    res = np.array(res)
    return res[:,0], res[:,1], res[:,2]
   
def compute_frac_and_save(filename, notes, wavelength, m_real, m_imag, r_min, r_max, nrad, a, df, rmon, nprocess=None):

    radii = np.logspace(np.log10(r_min),np.log10(r_max),nrad)
    wavelength_nm = wavelength*1e3

    assert nrad == 50
    nw = len(wavelength)

    w0_all = np.zeros((nrad,nw),np.float64)
    qext_all = np.zeros((nrad,nw),np.float64)
    g_all = np.zeros((nrad,nw),np.float64)

    for i,radius in enumerate(radii):
        print('Radius = %.4f um'%radius)
        qext, w0, g0 = fractal_meanfield(wavelength, m_imag, m_real, a, df, radius, rmon, nprocess)
        w0_all[i,:] = w0
        qext_all[i,:] = qext
        g_all[i,:] = g0

    # Interpolate to get rid of nans
    for i in range(len(radii)):
        inds = np.where(qext_all[i,:] != qext_all[i,:])[0]
        if len(inds) > 0:
            inds1 = np.where(qext_all[i,:] == qext_all[i,:])[0]
            qext_all[i,:] = np.interp(wavelength_nm,wavelength_nm[inds1],qext_all[i,inds1])
            w0_all[i,:] = np.interp(wavelength_nm,wavelength_nm[inds1],w0_all[i,inds1])
            g_all[i,:] = np.interp(wavelength_nm,wavelength_nm[inds1],g_all[i,inds1])

    utils.write_file(filename, notes, wavelength_nm, radii, w0_all, qext_all, g_all)