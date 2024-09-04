import utils
import numpy as np

def main():
    species = 'CO2'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    # Use Schmidt for longer wavelengths
    wv, xs = np.loadtxt('data/Schmidt2013_CO2.txt').T
    min_xs = np.min(wv)
    inds = np.where(out['wv']<np.min(wv))
    out['wv'] = np.append(out['wv'][inds],wv)
    out['xsa'] = np.append(out['xsa'][inds],xs)
    out['xsp'] = np.append(out['xsp'][inds],xs)
    out['xsi'] = np.append(out['xsi'][inds],xs*0.0)

    # Quantum yields
    ratios = {'wv': np.array([166.0, 168.0]),
        'CO2 + hv => CO + O': {'qy': np.array([0.0, 1.0])},
        'CO2 + hv => CO + O1D': {'qy': np.array([1.0, 0.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (phid, wogan,vulcan), ('Phidrates', 'Wogan','VULCAN'),xlim=(0,220))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Heays2017']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Schmidt2013']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Huebner2015']}
        ],
    'notes': "I use a low resolution version of the 270 K Schmidt et al. (2013) data."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()