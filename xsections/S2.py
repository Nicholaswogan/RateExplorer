import utils
import numpy as np

def main():
    species = 'S2'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    # Some adjustments
    for key in ['xsa','xsp','xsi']:
        for i in range(len(out[key])):
            if out[key][i] < 1e-45:
                out[key][i] = 0.0
    inds = np.where(out['xsa'] < out['xsp'] + out['xsi'])
    out['xsa'][inds] = out['xsp'][inds] + out['xsi'][inds]

    ratios = {
        'wv': np.array([np.min(out['wv']), np.max(out['wv'])]),
        'S2 + hv => S + S': {'qy': np.array([1.0, 1.0]),'new': False},
    }
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, wogan), ('VULCAN', 'Wogan'),xlim=(0,500))
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