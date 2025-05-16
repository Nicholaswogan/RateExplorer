import utils
import numpy as np

def main():
    species = 'NO3'
    out = utils.get_phidrates(species)

    vulcan = utils.get_VULCAN(species)

    for key in out:
        out[key] = out[key][:-1]

    min_xs = np.min(out['wv'])

    # Quantum yields
    ratios = {
        'wv': np.array([580., 581., 582., 583., 584., 585., 586., 587., 588., 589., 590.]),
        'NO3 + hv => NO2 + O': {'qy': np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])},
        'NO3 + hv => NO + O2': {'qy':  np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan,), ('VULCAN',), xlim=(0,800))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Huebner2015']},
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