import utils
import numpy as np

def get_yields():
    new_names = {
    'sCH2/H2': 'CH4 + hv => 1CH2 + H2',
    'CH3/H': 'CH4 + hv => CH3 + H',
    'CH2/H/H': 'CH4 + hv => CH2 + H + H',
    'CH/H2/H': 'CH4 + hv => CH + H2 + H',
    }
    qys = utils.get_phidrates_yields('CH4')
    qys_new = {'wv': qys['wv']}
    for key in qys:
        if key == 'wv':
            continue
        qys_new[new_names[key]] = qys[key]
    qys = qys_new
    gans = {
        'wv': np.array([118.2, 121.6]),
        'CH4 + hv => 1CH2 + H2': np.array([0.17,0.48]),
        'CH4 + hv => CH3 + H': np.array([0.26,0.42]),
        'CH4 + hv => CH2 + H + H': np.array([0.48,0.03]),
        'CH4 + hv => CH + H2 + H': np.array([0.09,0.07]),
    }
    inds = np.where(qys['wv'] < 97.7)
    qys_new = {}
    for key in qys:
        if key == 'wv':
            qys_new[key] = np.append(qys[key][inds],gans[key])
        else:
            qys_new[key] = {'qy': np.append(qys[key][inds],gans[key])}
    ratios = qys_new

    return ratios

def main():
    species = 'CH4'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    wv, xs = np.loadtxt('data/Lee2001_CH4.txt').T
    inds1 = np.where(out['wv'] < 142)
    inds2 = np.where(wv > 142)
    out['wv'] = np.append(out['wv'][inds1],wv[inds2])
    out['xsa'] = np.append(out['xsa'][inds1],xs[inds2])
    out['xsp'] = np.append(out['xsp'][inds1],xs[inds2])
    out['xsi'] = np.append(out['xsi'][inds1],xs[inds2]*0.0)

    wv_max = np.max(out['wv'])

    wv, xs = np.loadtxt('data/Karkoschka1994_CH4.txt').T
    wv = np.append(out['wv'][-1]*1.01,np.append(wv[0]*0.99,wv))
    xs = np.append(0.0,np.append(0.0,xs))
    out['wv'] = np.append(out['wv'], wv)
    out['xsa'] = np.append(out['xsa'],xs)
    out['xsp'] = np.append(out['xsp'],xs*0.0)
    out['xsi'] = np.append(out['xsi'],xs*0.0)

    # Quantum yields
    ratios = get_yields()
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, phid, wogan), ('VULCAN','Phidrates', 'Wogan'),xlim=(0,1000))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), 142.0], 'citations': ['Heays2017']},
        {'nm-range': [142.0, float(wv_max)], 'citations': ['Lee2001']},
        {'nm-range': [float(wv_max), float(np.max(out['wv']))], 'citations': ['Karkoschka1994']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), 97.7], 'citations': ['Heays2017']},
        {'nm-range': [97.7, float(np.max(out['wv']))], 'citations': ['Gans2011']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()