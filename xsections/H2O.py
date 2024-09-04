import utils
import numpy as np

def main():
    species = 'H2O'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
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

    # Tack on Ranjan
    wv, xs = np.loadtxt('data/Ranjan2020_H2O.txt').T
    inds1 = np.where(out['wv'] < 192.056)
    inds2 = np.where(wv > 192.056)
    out['wv'] = np.append(out['wv'][inds1],wv[inds2])
    out['xsa'] = np.append(out['xsa'][inds1],xs[inds2])
    out['xsp'] = np.append(out['xsp'][inds1],xs[inds2])
    out['xsi'] = np.append(out['xsi'][inds1],xs[inds2]*0.0)

    # Make sure absorption does not exceed photolysis and ionization
    inds = np.where(out['xsa'] < out['xsp'] + out['xsi'])
    out['xsa'][inds] = out['xsp'][inds] + out['xsi'][inds]

    # Quantum yields
    ratios = {
        'wv': np.array([105.0, 120.1, 120.2, 123.0, 123.1, 145.0, 145.1, 185.0]),
        'H2O + hv => OH + H': {'qy': np.array([0.89, 0.89, 0.78, 0.78, 0.89, 0.89, 1.0, 1.0])},
        'H2O + hv => H2 + O1D': {'qy': np.array([0.11, 0.11, 0.1, 0.1, 0.11, 0.11, 0.0, 0.0])},
        'H2O + hv => O + H + H': {'qy': np.array([0.0, 0.0, 0.12, 0.12, 0.0, 0.0, 0.0, 0.0])},
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
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), 192.056], 'citations': ['Heays2017']},
        {'nm-range': [192.056, float(np.max(out['wv']))], 'citations': ['Ranjan2020']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), 120.1], 'citations': ['Burkholder2020','Stief1975']},
        {'nm-range': [120.1, 123.1], 'citations': ['Burkholder2020','Slanger1982']},
        {'nm-range': [123.1, float(np.max(out['wv']))], 'citations': ['Burkholder2020','Stief1975']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()