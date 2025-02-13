import utils
from scipy import constants as const
import numpy as np
import yaml
from scipy import optimize
from photochem.utils._format import yaml, Loader, MyDumper, flowmap, FormatReactions_main
from matplotlib import pyplot as plt
import os

def objective(x, cp_298, cp_der_298, nr_dof_sp):
    coeffs = x
    cp_298_2 = utils.heat_capacity_shomate(coeffs, 298)
    cp_der_298_2 = utils.heat_capacity_shomate_derivative(coeffs, 298)
    
    # Values that match at 298 K
    res1 = np.array([
        (cp_298 - cp_298_2)/cp_298, # heat capacity matches
        (cp_der_298 - cp_der_298_2)/cp_der_298, # derivative matches (with normalization)
    ])
    
    # Next, we ensure the heat capacity relaxes to what it should be at really cold T
    TT = np.linspace(10,100,10)
    cps = np.array([utils.heat_capacity_shomate(coeffs, Ts) for Ts in TT])
    cp0 = (5/2)*const.R + (nr_dof_sp/2)*const.R
    res2 = (cps - cp0)/(cp0*len(TT))

    return np.append(res1,res2)

def compute_low_T_thermo(species, nr_dof):

    sp = species['name']
    Tr = species['thermo']['temperature-ranges']

    assert Tr[0] <= 298
    assert Tr[1] >= 298

    coeffs1 = species['thermo']['data'][0]
    cp_298 = utils.heat_capacity_shomate(coeffs1, 298)
    cp_der_298 = utils.heat_capacity_shomate_derivative(coeffs1, 298)
    nr_dof_sp = nr_dof[sp]

    guess = [cp_298, 0.0, 0.0, 0.0, 0.0]
    args = (cp_298, cp_der_298, nr_dof_sp,)
    sol = optimize.root(objective, guess, args=args, method='lm')
    if not sol.success:
        print(sp)

    coeffs = sol.x

    H_298 = utils.enthalpy_shomate(coeffs1, 298)/1000
    S_298 = utils.entropy_shomate(coeffs1, 298)

    TT = 298/1000
    a5 = H_298 - (coeffs[0]*TT + (coeffs[1]*TT**2)/2.0 \
                + (coeffs[2]*TT**3)/3.0  + (coeffs[3]*TT**4)/4.0 \
                - coeffs[4]/TT)

    a6 = S_298 - (coeffs[0]*np.log(TT) + coeffs[1]*TT \
                + (coeffs[2]*TT**2)/2.0 + (coeffs[3]*TT**3)/3.0 \
                - coeffs[4]/(2.0 * TT**2))

    coeffs = np.append(coeffs,[a5,a6])

    return coeffs

def main():
    with open('thermodata121.yaml','r') as f:
        mech = yaml.load(f,yaml.Loader)
    with open('rotational_dof.yaml','r') as f:
        nr_dof = yaml.load(f,yaml.Loader)

    new_data = {}
    for i,species in enumerate(mech['species']):
        if species['name'] in ['He','H']:
            continue
        coeffs = compute_low_T_thermo(species, nr_dof)
        coeffs = [float(a) for a in coeffs]
        new_data[species['name']] = coeffs

    with open('tmp.yaml','w') as f:
        yaml.dump(new_data,f,sort_keys=False,Dumper=utils.CustomDumper,width=75)
    with open('tmp.yaml','r') as f:
        new_data = yaml.load(f,yaml.Loader)
    os.remove('tmp.yaml')

    for i,species in enumerate(mech['species']):
        if species['name'] not in new_data:
            continue
        new_Tr = [10.0, 298.0]+ mech['species'][i]['thermo']['temperature-ranges'][1:]
        new_coeffs = [new_data[species['name']]] + mech['species'][i]['thermo']['data']

        mech['species'][i]['thermo']['temperature-ranges'] = new_Tr
        mech['species'][i]['thermo']['data'] = new_coeffs

    mech = FormatReactions_main(mech)
    with open('thermodata122.yaml','w') as f:
        yaml.dump(mech,f,sort_keys=False,Dumper=yaml.Dumper,width=75)

def plot(sp, nr_dof, gas):
    fig,axs = plt.subplots(1,3,figsize=[11,3])

    ax = axs[0]

    TT = np.linspace(298,350,100)
    ax.plot(TT, [gas.heat_capacity(sp, T) for T in TT])

    TT = np.linspace(10,298,100)
    ax.plot(TT, [gas.heat_capacity(sp, T) for T in TT])

    cp0 = (5/2)*const.R + (nr_dof[sp]/2)*const.R
    ax.axhline(cp0, c='k', ls='--')
    ax.set_ylabel('Heat capacity (J/mol/K)')
    ax.set_xlabel('Temperature (K)')


    ax = axs[1]

    TT = np.linspace(298,310,100)
    ax.plot(TT, [gas.enthalpy(sp, T)/1e3 for T in TT])

    TT = np.linspace(280,298,100)
    ax.plot(TT, [gas.enthalpy(sp, T)/1e3 for T in TT])

    ax.set_ylabel('Enthalpy (kJ/mol)')
    ax.set_xlabel('Temperature (K)')

    ax = axs[2]

    TT = np.linspace(298,310,100)
    ax.plot(TT, [gas.entropy(sp, T) for T in TT])

    TT = np.linspace(280,298,100)
    ax.plot(TT, [gas.entropy(sp, T) for T in TT])

    ax.set_ylabel('Entropy (J/mol)')
    ax.set_xlabel('Temperature (K)')

    for ax in axs:
        ax.grid(alpha=0.4)

    plt.subplots_adjust(wspace=0.4)

    plt.savefig('figures/'+sp+'.pdf',bbox_inches = 'tight')

def plot_all():
    with open('rotational_dof.yaml','r') as f:
        nr_dof = yaml.load(f,yaml.Loader)
    gas = utils.ReactionExplorer('thermodata122.yaml')

    for sp in gas.mech['species']:
        plot(sp, nr_dof, gas)


if __name__ == '__main__':
    main()
    # plot_all()


