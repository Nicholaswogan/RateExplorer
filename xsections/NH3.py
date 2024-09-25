import utils
import numpy as np

def main():
    species = 'NH3'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    # Prepend phidrates
    min_xs = np.min(out['wv'])
    out = utils.prepend_xs_of_other(out, phid)

    # Append Cheng2006
    max_xs = np.max(out['wv'])
    wv, xs = np.loadtxt('data/Cheng2006_NH3.txt').T
    inds = np.where(wv > out['wv'][-1])
    out['wv'] = np.append(out['wv'], wv[inds])
    out['xsp'] = np.append(out['xsp'], xs[inds])
    out['xsa'] = np.append(out['xsa'], xs[inds])
    out['xsi'] = np.append(out['xsi'], xs[inds]*0)

    # Quantum yields
    new_names = {
        'NH2/H': 'NH3 + hv => NH2 + H',
        'sNH/H2': 'NH3 + hv => NH + H2',
        'NH/H/H': 'NH3 + hv => NH + H + H',
    }
    qys = utils.get_phidrates_yields('NH3')
    ratios = {'wv': qys['wv']}
    _, inds_unique = np.unique(ratios['wv'].astype(np.float32),return_index=True) # Remove duplicate entries
    ratios['wv'] = ratios['wv'][inds_unique]
    for key in qys:
        if key == 'wv':
            continue
        ratios[new_names[key]] = {'qy': qys[key][inds_unique]}
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, phid, wogan), ('VULCAN','Phidrates', 'Wogan'),xlim=(0,300))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(max_xs)], 'citations': ['Heays2017']},
        {'nm-range': [float(max_xs), float(np.max(out['wv']))], 'citations': ['KellerRudek2013','Cheng2006']}
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