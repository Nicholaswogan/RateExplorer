import utils
import numpy as np

def main():
    species = 'HS'
    out = utils.get_leiden('SH')

    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    # Make sure absorption does not exceed photolysis and ionization
    inds = np.where(out['xsa'] < out['xsp'] + out['xsi'])
    out['xsa'][inds] = out['xsp'][inds] + out['xsi'][inds]

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'HS + hv => S + H': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan,wogan,), ('VULCAN','Wogan',))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()