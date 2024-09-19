import utils
import numpy as np

def main():
    species = 'C3H6'

    wv, xs = np.loadtxt('data/Walker2007_C3H6.txt').T
    out = {}
    out['wv'] = wv
    out['xsa'] = xs
    out['xsp'] = xs
    out['xsi'] = xs*0.0

    # Quantum yields
    ratios = {
        'wv': np.array([100.0, 135.0, 136.0, 155.0, 156.0, 175.0, 176.0, 195.0]),
        'C3H6 + hv => C3H4 + H2': {'qy': np.array([0.43+0.25, 0.43+0.25, 0.40+0.24, 0.40+0.24, 
                                                   0.015+0.015+0.565, 0.015+0.015+0.565,
                                                   0.01+0.01+0.6, 0.01+0.01+0.6])},
        'C3H6 + hv => C2H4 + 1CH2': {'qy': np.array([0.06, 0.06, 0.04, 0.04, 0.02, 0.02, 0.0, 0.0])},
        'C3H6 + hv => C2H3 + CH3': {'qy': np.array([0.21, 0.21, 0.27, 0.27, 0.335, 0.335, 0.34, 0.34])},
        'C3H6 + hv => C2H2 + CH4': {'qy': np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, xlim=(0,200))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Walker2007']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Lavvas2008']}
        ],
    'notes': "The Lavvas et al. (2008) paths making CH2CCH2 + H2, CH3C2H + H2 and C3H5 + H are all in the C3H4 + H2 branch."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()