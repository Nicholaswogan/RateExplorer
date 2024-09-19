import utils
import numpy as np

def main():
    species = 'OClO'
    out = utils.get_phidrates(species)

    # Quantum yields
    ratios = {
        'wv': np.array([200.0, 276.0, 277.0, 400.0]),
        'OClO + hv => ClO + O1D': {'qy': np.array([1.0, 1.0, 0.0, 0.0])},
        'OClO + hv => ClO + O': {'qy': np.array([0.0, 0.0, 1.0, 1.0])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, xlim=(0,600))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Huebner2015']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Huebner2015']}
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()