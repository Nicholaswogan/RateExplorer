import utils
import numpy as np

def main():
    species = 'C4H4'
    wogan = utils.get_wogan(species)

    wv, xs = np.loadtxt('data/Fahr1996_C4H4.txt').T
    out = {}
    out['wv'] = wv
    out['xsa'] = xs
    out['xsp'] = xs
    out['xsi'] = xs*0

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'C4H4 + hv => C4H2 + H2': {'qy': np.array([0.8, 0.8])},
        'C4H4 + hv => C2H2 + C2H2': {'qy': np.array([0.2, 0.2])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,), ('Wogan',),xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Fahr1996']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Lavvas2008']}
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()