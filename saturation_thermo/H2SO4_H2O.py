import utils
import numpy as np
import ctypes as ct
import numba as nb
import subprocess
from matplotlib import pyplot as plt
from zeleznik import saturation_vapor_pressures

def create_binary_saturation_pressure():
    subprocess.call('gfortran sulfuric_acid.f90 -fPIC -O3 -o sulfuric_acid.so'.split())
    binary_saturation_pressure_c = ct.CDLL('sulfuric_acid.so').binary_saturation_pressure_c
    binary_saturation_pressure_c.argtypes = [ct.c_double, ct.c_double, ct.c_void_p, ct.c_void_p]
    binary_saturation_pressure_c.restype = None
    
    @nb.njit(nb.types.Array(nb.float64,1,'C')(nb.float64,nb.float64))
    def binary_saturation_pressure(T, x_H2SO4):
        P_H2SO4, P_H2O = np.array(0.0,np.double), np.array(0.0,np.double)
        binary_saturation_pressure_c(T, x_H2SO4, P_H2SO4.ctypes.data, P_H2O.ctypes.data)
        return np.array([P_H2SO4.item(), P_H2O.item()])
        
    return binary_saturation_pressure

binary_saturation_pressure = create_binary_saturation_pressure()

def print_variable(arr, varname):
    str = 'real(dp), parameter :: '+varname+'(*) = &\n'
    str += '[ &\n'
    i = 0
    counter = 0
    while True:
    
        if i >= len(arr):
            str = str[:-2]
            str += ' &'
            break
    
        # print('%e, '%(arr[i]),end='')
        str += '%e_dp, '%(arr[i])
        
        counter += 1
        if counter > 4:
            counter = 0
            # print('&')
            str += '&\n'
    
        i += 1
    str += '\n]'
    print(str)

def fit_to_zeleznik():
    h_H2O = utils.SaturationProperties()
    h_H2O.init2('results/H2O_sat.yaml')

    h_H2SO4 = utils.SaturationProperties()
    h_H2SO4.init2('results/H2SO4_sat.yaml')

    x1s = np.arange(0.05,1.0-1e-10,0.05)
    x1s = np.append([1e-4,0.01],x1s)
    x1s = np.append(x1s,[0.99,1-1e-4])

    res_H2SO4 = {}
    res_H2O = {}
    keys = ['P_ref','A_v','B_v','A_s','B_s','A_c','B_c']
    for key in keys:
        res_H2SO4[key] = np.zeros(x1s.shape[0])
        res_H2O[key] = np.zeros(x1s.shape[0])

    for i,x1 in enumerate(x1s):
        mu = 18.01534
        T_triple = 273.15
        T_critical = 647.0
        T_ref = 300.0 # K
        P_ref_H2SO4, P_ref_H2O = saturation_vapor_pressures(T_ref, x1, h_H2SO4, h_H2O)
        P_ref_H2SO4 = P_ref_H2SO4.item()
        P_ref_H2O = P_ref_H2O.item()
        
        T_data = np.linspace(200.0, 350.0, 20)
        T_data_H2SO4 = T_data.copy()
        T_data_H2O = T_data.copy()
        P_data = np.array([saturation_vapor_pressures(T, x1, h_H2SO4, h_H2O) for T in T_data])
        P_data_H2SO4 = P_data[:,0].astype(np.float64)
        P_data_H2O = P_data[:,1].astype(np.float64)

        if i == 0:
            P_ref_H2O = h_H2O.sat_pressure(T_ref)
            T_data_H2O = np.linspace(100.0, 600.0, 30)
            P_data_H2O = np.array([h_H2O.sat_pressure(T) for T in T_data_H2O])

        if i == x1s.shape[0]-1:
            P_ref_H2SO4 = h_H2SO4.sat_pressure(T_ref)
            T_data_H2SO4 = np.linspace(200.0, 600.0, 30)
            P_data_H2SO4 = np.array([h_H2SO4.sat_pressure(T) for T in T_data_H2SO4])
        
        h1 = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref_H2SO4)
        guess = [4e9, 0, 6e9, 0]
        sol = h1.optimize(guess, T_data_H2SO4, P_data_H2SO4)
        
        h2 = utils.SaturationPropertiesFitter(mu, T_triple, T_critical, T_ref, P_ref_H2O)
        guess = [4e9, 0, 6e9, 0]
        sol = h2.optimize(guess, T_data_H2O, P_data_H2O)

        res_H2SO4['P_ref'][i] = h1.P_ref
        res_H2SO4['A_v'][i] = h1.A_v
        res_H2SO4['B_v'][i] = h1.B_v
        res_H2SO4['A_s'][i] = h1.A_s
        res_H2SO4['B_s'][i] = h1.B_s
        res_H2SO4['A_c'][i] = h1.A_c
        res_H2SO4['B_c'][i] = h1.B_c

        res_H2O['P_ref'][i] = h2.P_ref
        res_H2O['A_v'][i] = h2.A_v
        res_H2O['B_v'][i] = h2.B_v
        res_H2O['A_s'][i] = h2.A_s
        res_H2O['B_s'][i] = h2.B_s
        res_H2O['A_c'][i] = h2.A_c
        res_H2O['B_c'][i] = h2.B_c

    out = {}
    out['mu'] = 18.01534
    out['T-triple'] = 273.15
    out['T-critical'] = 647.0
    out['T-ref'] = 300.0
    out['x_H2SO4'] = x1s
    out['x_H2SO4'][0] = 0.0
    out['x_H2SO4'][-1] = 1
    out['H2SO4'] = res_H2SO4
    out['H2O'] = res_H2O

    arr = out['x_H2SO4']
    varname = 'x_H2SO4_grid'
    print_variable(arr, varname)
    print()
    tmp = 'H2SO4'
    arr = out[tmp]['P_ref']
    varname = 'P_ref_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_v']
    varname = 'A_v_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_v']
    varname = 'B_v_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_s']
    varname = 'A_s_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_s']
    varname = 'B_s_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_c']
    varname = 'A_c_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_c']
    varname = 'B_c_'+tmp
    print_variable(arr, varname)
    print()
    tmp = 'H2O'
    arr = out[tmp]['P_ref']
    varname = 'P_ref_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_v']
    varname = 'A_v_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_v']
    varname = 'B_v_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_s']
    varname = 'A_s_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_s']
    varname = 'B_s_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['A_c']
    varname = 'A_c_'+tmp
    print_variable(arr, varname)
    arr = out[tmp]['B_c']
    varname = 'B_c_'+tmp
    print_variable(arr, varname)

def plot_comparison():
    h_H2O = utils.SaturationProperties()
    h_H2O.init2('results/H2O_sat.yaml')

    h_H2SO4 = utils.SaturationProperties()
    h_H2SO4.init2('results/H2SO4_sat.yaml')

    plt.rcParams.update({'font.size': 12})
    fig,axs = plt.subplots(1,2,figsize=[10,5])

    ax = axs[0]
    ax1 = axs[1]
    TT = np.linspace(150,350,15)
    TT1 = np.linspace(150,600,30)
    x_H2SO4 = np.array([1.0e-4, 0.075, 0.511, 0.98])

    ax.plot(TT1,[h_H2SO4.sat_pressure(T)/1e6 for T in TT1], c='k', lw=2)
    ax1.plot(TT1,[h_H2O.sat_pressure(T)/1e6 for T in TT1], c='k', lw=2)
    ax.plot([],[],c='k',ls='-', lw=2, label='Pure compositions')

    lss = ['--',':','-.','']
    markers= ['','','','*']
    for i in range(x_H2SO4.shape[0]):
        x1 = x_H2SO4[i]
        ls = lss[i]
        marker = markers[i]
        tmp = np.array([saturation_vapor_pressures(T, x1, h_H2SO4, h_H2O)/1e6 for T in TT])
        ax.plot(TT,tmp[:,0], c='C0', ls=ls, lw=2,marker=marker)
        ax1.plot(TT,tmp[:,1], c='C0', ls=ls, lw=2,marker=marker)
        tmp1 = np.array([binary_saturation_pressure(T, x1) for T in TT1])
        ax.plot(TT1,tmp1[:,0]/1e6, c='C1', ls=ls, lw=1,marker=marker)
        ax1.plot(TT1,tmp1[:,1]/1e6, c='C1', ls=ls, lw=1,marker=marker)
        ax.plot([],[],c='k',ls=ls, lw=2, label='$x_\mathrm{H_2SO_4}$ = %.2e'%(x1),marker=marker)

    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.plot([],[],c='C0',lw=2,label='Zeleznik (1991)')
    ax2.plot([],[],c='C1',lw=2,label='sulfuric_acid.f90 (interpolation)')

    # ax1.legend()
    ax.legend(ncol=1,bbox_to_anchor=(0, 1.02), loc='lower left')
    ax2.legend(ncol=1,bbox_to_anchor=(1.6, 1.02), loc='lower right')

    ax.set_yscale('log')
    ax.set_ylim(1e-30,1)
    ax.set_ylabel('Saturation vapor pressure (bar)')
    ax.set_xlabel('Temperature (K)')
    ax.grid(alpha=0.4)
    ax.text(0.02, .98, 'H$_2$SO$_4$', \
            size = 20,ha='left', va='top',transform=ax.transAxes)

    ax1.set_yscale('log')
    ax1.set_ylim(1e-25,3e2)
    # ax1.set_ylabel('Saturation vapor pressure (bar)')
    ax1.set_xlabel('Temperature (K)')
    ax1.grid(alpha=0.4)
    ax1.text(0.02, .98, 'H$_2$O', \
            size = 20,ha='left', va='top',transform=ax1.transAxes)

    plt.savefig('figures/Zeleznik_1991_vs_interpolation.pdf',bbox_inches='tight')


def main():
    fit_to_zeleznik()
    plot_comparison()

if __name__ == '__main__':
    main()