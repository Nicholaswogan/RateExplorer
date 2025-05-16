import utils
import numpy as np

def main():
    species = 'CN'
    out = utils.get_leiden(species)
    phi = utils.get_phidrates(species)

    min_xs = np.min(out['wv'])

    inds = np.where(phi['wv']<out['wv'][0])
    out['wv'] = np.append(phi['wv'][inds],out['wv'])
    out['xsa'] = np.append(phi['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(phi['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(phi['xsa'][inds]*0,out['xsp'])

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'CN + hv => C + N': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (phi,), ('Phidrates',), xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']},
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