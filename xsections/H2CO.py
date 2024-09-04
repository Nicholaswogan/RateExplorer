import utils
import numpy as np

def get_ratios():

    # Get Phidrates quantum yields
    phi = utils.phidrates('H2CO')
    xsp = np.zeros(phi.neutral.data['wavelength'].shape)
    for branch in phi.neutral.branches:
        xsp += phi.neutral.data[branch]
    xsp = phi.neutral.data['cross section']*xsp
    qys = {'wv':phi.neutral.data['wavelength']}
    for branch in phi.neutral.branches:
        qy = (phi.neutral.data[branch]*phi.neutral.data['cross section'])/xsp
        qys[branch] = qy

    # Compute the qy of H + HCO based on JPL-19
    wv = qys['wv']
    a0 = 557.95835182
    a1 = -7.31994058026
    a2 = 0.03553521598
    a3 = -7.54849718e-5
    a4 = 5.91001021e-8
    phi_1 = a0 + a1*wv + a2*wv**2 + a3*wv**3 + a4*wv**4
    phi_1[qys['wv'] > 338] = 0.0
    phi_1[qys['wv'] < 255] = 0.0
    phi_1[(qys['wv'] < 255) & (qys['wv'] > 198)] = 0.3

    # Assume the Phidrates quantum yields fill in the rest
    phi_tot = 1 - phi_1
    phi_2 = phi_tot*qys['CO/H2']/(qys['CO/H/H']+qys['CO/H2'])
    phi_3 = phi_tot*qys['CO/H/H']/(qys['CO/H/H']+qys['CO/H2'])

    ratios = {
        'wv': wv,
        'H2CO + hv => HCO + H': {'qy': phi_1, 'new': False},
        'H2CO + hv => CO + H2': {'qy': phi_2, 'new': False},
        'H2CO + hv => CO + H + H': {'qy': phi_3, 'new': False}
    }

    return ratios

def main():
    species = 'H2CO'
    out = utils.get_leiden(species)
    ld = utils.get_leiden(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    min_xs = np.min(out['wv'])
    
    # Prepend Phidrates assuming its all ionization
    other = phid
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds]*0,out['xsp'])

    # Get ratios
    out['ratios'] = get_ratios()
    out['missing'] = []

    # Make plots
    utils.make_xs_plot(species, out, (ld, phid), ('Leiden','Phidrates'),xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), 250.0], 'citations': ['Huebner2015']},
        {'nm-range': [250.0, float(np.max(out['wv']))], 'citations': ['Huebner2015','Burkholder2020']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()