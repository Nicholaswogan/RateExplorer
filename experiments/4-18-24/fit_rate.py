import numpy as np
from matplotlib import pyplot as plt
import utils
from scipy import optimize
import pandas as pd

def objective(x, T_data, k_data):
    A, b, Ea = x
    k_ = utils.arrhenius_rate(A, b, Ea, T_data)
    return np.log10(k_) - np.log10(k_data)

def fit_rate_constant():
    r = utils.ReactionExplorer('thermodata121_test.yaml')
    P = 1 # bars
    T_data1 = np.array([190,400,30])
    k_data1 = np.empty(T_data1.shape[0])
    for i in range(T_data1.shape[0]):
        tmp  = r.evaluate_rates(T_data1[i], P)
        k_data1[i] = tmp['kircher1']['k_f'] + tmp['kircher2']['k_f']
    
    inds = np.where(T_data1<200)
    k_data1[inds] = 1.3e-11
    
    T_data2 = np.array([401,1500,30])
    k_data2 = np.empty(T_data2.shape[0])
    for i in range(T_data2.shape[0]):
        tmp  = r.evaluate_rates(T_data2[i], P)
        k_data2[i] = tmp['klippenstein1']['k_f']
    
    T_data = np.append(T_data1,T_data2)
    k_data = np.append(k_data1,k_data2)

    sol = optimize.root(objective, [2.83e-20, 2, -1304], method='lm', args=(T_data,k_data))
    
    A, b, Ea = sol.x
    print('- equation: HO2 + HO2 <=> H2O2 + O2')
    print('  rate-constant: {A: %.2e, b: %.2f, Ea: %.1f}'%(A,b,Ea))
    return A, b, Ea

def plot(A, b, Ea):
    r = utils.ReactionExplorer('thermodata121_test.yaml')
    res = pd.read_excel('data/ktp_wogan.xlsx')

    plt.rcParams.update({'font.size': 13})
    fig,axs = plt.subplots(1,2,figsize=[10,4],sharey=True)

    ax = axs[0]
    T = np.linspace(150,400,100)
    P = 1
    k1 = np.empty(T.shape[0])
    k2 = np.empty(T.shape[0])
    for i,T1 in enumerate(T):
        tmp  = r.evaluate_rates(T1, P)
        k1[i] = utils.arrhenius_rate(A, b, Ea, T1)
        k2[i] = tmp['wogan']['k_f']
    assert np.all(np.isclose(k1,k2))
    ax.plot(T, k1, c='C2', ls='-', label='%.2e * T^%.2f * exp(-(%.1f)/T)'%(A, b, Ea))
    # ax.plot(T, k2, c='C2', ls='-', label='%.2e * T^%.2f * exp(-(%.1f)/T)'%(A, b, Ea))

    P = np.array([1, 0.1, 0])
    ls = ['-',':','--']
    for j,P1 in enumerate(P):
        k1 = np.empty(T.shape[0])
        for i,T1 in enumerate(T):
            tmp  = r.evaluate_rates(T1, P1)
            k1[i] = tmp['kircher1']['k_f'] + tmp['kircher2']['k_f']
        ax.plot(T, k1, c='C0', ls=ls[j], label='Kircher %.1f bar'%(P1))

    P = 0
    k1 = np.empty(T.shape[0])
    for i,T1 in enumerate(T):
        tmp  = r.evaluate_rates(T1, P)
        k1[i] = tmp['klippenstein1']['k_f']
    ax.plot(T, k1, c='C1', ls='-', label='Klippenstein fit (1 bar)')

    ms = 3
    Tk = [res[a][44] for a in res][1:]
    ratek = [res[a][52] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='-',marker='o', ms=ms, label='Klippenstein 0.1 bar')

    ratek = [res[a][54] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='--',marker='o', ms=ms, label='Klippenstein 1 bar')

    ratek = [res[a][56] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='-.',marker='o', ms=ms, label='Klippenstein 10 bar')

    ax.set_xlim(150,400)
    ax.set_yscale('log')
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Rate constant (cm$^3$ molecules$^{-1}$ s$^{-1}$)')
    ax.grid(alpha=0.4)
    ax.legend(ncol=3,bbox_to_anchor=(0, 1.02), loc='lower left',fontsize=11)

    ax = axs[1]
    T = np.linspace(150,2000,100)

    P = 1
    k1 = np.empty(T.shape[0])
    k2 = np.empty(T.shape[0])
    for i,T1 in enumerate(T):
        tmp  = r.evaluate_rates(T1, P)
        k1[i] = utils.arrhenius_rate(A, b, Ea, T1)
        k2[i] = tmp['wogan']['k_f']
    assert np.all(np.isclose(k1,k2))
    ax.plot(T, k1, c='C2', ls='-')
    # ax.plot(T, k2, c='C2', ls='-')

    P = np.array([1, 0.1, 0])
    ls = ['-',':','--']
    for j,P1 in enumerate(P):
        k1 = np.empty(T.shape[0])
        for i,T1 in enumerate(T):
            tmp  = r.evaluate_rates(T1, P1)
            k1[i] = tmp['kircher1']['k_f'] + tmp['kircher2']['k_f']
        ax.plot(T, k1, c='C0', ls=ls[j], label='Kircher %.1f bar'%(P1))

    P = 0
    k1 = np.empty(T.shape[0])
    for i,T1 in enumerate(T):
        tmp  = r.evaluate_rates(T1, P)
        k1[i] = tmp['klippenstein1']['k_f']
    ax.plot(T, k1, c='C1', ls='-', label='Klippenstein')

    Tk = [res[a][44] for a in res][1:]
    ratek = [res[a][52] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='-',marker='o', ms=ms, label='Klippenstein 0.1 bar')

    ratek = [res[a][54] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='--',marker='o', ms=ms, label='Klippenstein 1 bar')

    ratek = [res[a][56] for a in res][1:]
    ax.plot(Tk, ratek, c='k', ls='-.',marker='o', ms=ms, label='Klippenstein 10 bar')

    ax.set_yscale('log')
    ax.set_xlabel('Temperature (K)')
    ax.grid(alpha=0.4)

    plt.subplots_adjust(wspace=0.05)
    plt.savefig('new_rate.pdf',bbox_inches='tight')

if __name__ == '__main__':
    A, b, Ea = fit_rate_constant()
    plot(A, b, Ea)