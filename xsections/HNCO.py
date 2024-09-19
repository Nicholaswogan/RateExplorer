import utils
import numpy as np

def main():
    species = 'HNCO'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phi = utils.get_phidrates(species)

    # Prepend Phidrates assuming its all ionization
    min_xs = np.min(out['wv'])
    other = phi
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds]*0,out['xsp'])

    # Quantum yields
    ratios = {
        'wv': np.array([100.0, 250, 255.0, 300.0]),
        'HNCO + hv => NH + CO': {'qy': np.array([0.5, 0.5, 1.0, 1.0])},
        'HNCO + hv => H + NCO': {'qy': np.array([0.5, 0.5, 0.0, 0.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan,phi), ('VULCAN','Phidrates'), xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']},
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