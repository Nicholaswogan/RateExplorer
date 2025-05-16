import utils
import numpy as np

def main():
    species = 'ClO'
    out = utils.get_wogan(species)
    del out['ratios']
    del out['missing']

    min_xs = np.min(out['wv'])

    _, inds_unique = np.unique(out['wv'].astype(np.float32),return_index=True)
    for key in out:
        out[key] = out[key][inds_unique]

    # Quantum yields
    ratios = {
        'wv': np.array([263.3, 263.5]),
        'ClO + hv => Cl + O1D': {'qy': np.array([1.0, 0.0])},
        'ClO + hv => Cl + O': {'qy': np.array([0.0, 1.0])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (), (), xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Schmidt1998']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Burkholder1990']}
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()