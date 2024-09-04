import utils
import numpy as np

def main():
    species = 'NO2'
    out = utils.get_phidrates(species)
    
    ld = utils.get_leiden(species)
    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    with open('data/JPL19_NO2.txt') as f:
        lines = f.readlines()
    line = lines[0]
    wv = (np.array([float(a.split('–')[0]) for a in line.split()]) + np.array([float(a.split('–')[1]) for a in line.split()]))/2
    xs = np.array([float(a) for a in lines[1].split()])*1e-20

    max_wv = out['wv'][-1]
    inds = np.where(wv > out['wv'][-1])
    out['wv'] = np.append(out['wv'],wv[inds])
    out['xsa'] = np.append(out['xsa'],xs[inds])
    out['xsp'] = np.append(out['xsp'],xs[inds]*0)
    out['xsi'] = np.append(out['xsi'],xs[inds]*0)

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'NO2 + hv => NO + O': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, wogan, ld), ('VULCAN', 'Wogan', 'Leiden'),xlim=(0,700))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(max_wv)], 'citations': ['Huebner2015']},
        {'nm-range': [float(max_wv), float(np.max(out['wv']))], 'citations': ['Burkholder2020']},
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