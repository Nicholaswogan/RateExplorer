import utils
import numpy as np

def main():
    species = 'CH3'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)

    for key in ['xsa','xsp','xsi']:
        inds = np.where(out[key] < 1e-45)
        out[key][inds] = 0.0

    # Make sure absorption does not exceed photolysis and ionization
    inds = np.where(out['xsa'] < out['xsp'] + out['xsi'])
    out['xsa'][inds] = out['xsp'][inds] + out['xsi'][inds]

    # Quantum yields
    ratios = {
        'wv': np.array([164.9, 165.1]),
        'CH3 + hv => 1CH2 + H': {'qy': np.array([1.0, 0.0])},
        'CH3 + hv => CH2 + H': {'qy': np.array([0.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan,), ('VULCAN',), xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']},
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