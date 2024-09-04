import utils
import numpy as np

def main():
    species = 'C2H6'
    out = utils.get_leiden(species)

    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)
    
    # Prepend Phidrates assuming its all ionization
    min_xs = np.min(out['wv'])
    other = phid
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds]*0,out['xsp'])

    # Yields
    ratios = {
        'wv': np.array([120.6, 121.6, 122.6]),
        'C2H6 + hv => CH3 + CH3': {'qy': np.array([0.17, 0.03, 0.0])},
        'C2H6 + hv => CH4 + 1CH2': {'qy': np.array([0.16, 0.26, 0.02])},
        'C2H6 + hv => C2H4 + H + H': {'qy': np.array([0.41, 0.31, 0.13])},
        'C2H6 + hv => C2H4 + H2': {'qy': np.array([0.0, 0.15, 0.48])},
        'C2H6 + hv => C2H2 + H2 + H2': {'qy': np.array([0.26, 0.25, 0.37])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (phid, wogan), ('Phidrates', 'Wogan'),xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']}
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