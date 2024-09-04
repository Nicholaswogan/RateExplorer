import utils
import numpy as np

def main():
    species = 'SO2'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)    

    # Quantum yields
    ratios = {
        'wv': np.array([206.5, 207.0]),
        'SO2 + hv => SO + O': {'qy': np.array([0.5, 1.0])},
        'SO2 + hv => S + O2': {'qy': np.array([0.5, 0.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, phid, wogan), ('VULCAN','Phidrates', 'Wogan'),xlim=(0,500))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
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