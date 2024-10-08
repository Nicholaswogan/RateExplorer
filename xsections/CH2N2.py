import utils
import numpy as np

def main():
    species = 'CH2N2'
    wv, xs = np.loadtxt('data/McMillan1966_CH2N2.txt').T
    out = {}
    out['wv'] = wv
    out['xsp'] = xs
    out['xsa'] = xs
    out['xsi'] = xs*0   

    wogan = utils.get_wogan(species)

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'CH2N2 + hv => 1CH2 + N2': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,), ('Wogan',),xlim=None)
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['KellerRudek2013']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Assumed']}
        ],
    'notes': 'I used the T = 298 K file.'
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()