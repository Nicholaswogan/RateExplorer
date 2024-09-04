import utils
import numpy as np

def get_xs():
    wv1, xs1 = np.loadtxt('data/BiehlStuhl1991_N2H4.txt').T
    wv2, xs2 = np.loadtxt('data/Vaghjiani1993_N2H4.txt').T
    inds = np.where(wv2 > wv1[-1])
    wv = np.append(wv1,wv2[inds])
    xs = np.append(xs1,xs2[inds])
    out = {}
    out['wv'] = wv
    out['xsp'] = xs
    out['xsa'] = xs
    out['xsi'] = xs*0
    return out

def main():
    species = 'N2H4'
    out = get_xs()

    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    # Quantum yields
    ratios = {
        'wv': utils.minmaxarray(out['wv']),
        'N2H4 + hv => N2H3 + H': {'qy': np.array([1.0, 1.0])}
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan,wogan,), ('VULCAN','Wogan',),xlim=None, ylim=None)
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), 239.945], 'citations': ['KellerRudek2013','BiehlStuhl1991']},
        {'nm-range': [239.945, float(np.max(out['wv']))], 'citations': ['KellerRudek2013','Vaghjiani1993']}
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