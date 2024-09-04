import utils
import numpy as np

def main():
    species = 'CO'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates('sCO')
    wogan = utils.get_wogan(species)

    # Where photolysis xs beocmes non-zero, prepend Phidrates
    for i in range(len(out['wv'])):
        if out['xsp'][i] > 0:
            ind = i
            min_xs = out['wv'][i]
            break

    other = phid
    inds = np.where(other['wv']<out['wv'][ind])
    out['wv'] = np.append(other['wv'][inds],out['wv'][ind:])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'][ind:])
    out['xsi'] = np.append(other['xsi'][inds],out['xsi'][ind:])
    out['xsp'] = np.append(other['xsp'][inds],out['xsp'][ind:])

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'CO + hv => C + O': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)
    
    # Make plots
    utils.make_xs_plot(species, out, (phid, wogan, vulcan), ('Phidrates', 'Wogan','VULCAN'),xlim=(0,200))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']}
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