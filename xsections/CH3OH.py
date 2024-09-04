import utils
import numpy as np

def main():
    species = 'CH3OH'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'CH3OH + hv => CH3O + H': {'qy':  np.array([0.86, 0.86])},
        'CH3OH + hv => CH3 + OH': {'qy':  np.array([0.07, 0.07])},
        'CH3OH + hv => H2CO + H2': {'qy':  np.array([0.07, 0.07])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, phid, wogan), ('VULCAN','Phidrates', 'Wogan'),xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Hagege1968']}
        ],
    'notes': "Hagege et al. reports two quantum yields, and I average them."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()