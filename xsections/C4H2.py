import utils
import numpy as np

def main():
    species = 'C4H2'

    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    wv1,xs1 = np.loadtxt('data/Ferradaz2009_C4H2.txt').T
    wv2,xs2 = np.loadtxt('data/Smith1998_C4H2.txt').T
    inds = np.where(wv1 < wv2[0])
    wv = np.append(wv1[inds],wv2)
    xs = np.append(xs1[inds],xs2)
    out = {}
    out['wv'] = wv
    out['xsp'] = xs
    out['xsa'] = xs
    out['xsi'] = xs*0.0

    # Quantum yields
    ratios = {
        'wv': np.array([140.0, 164.0, 165.0, 204.0, 205.0, 240.0]),
        'C4H2 + hv => C4H + H': {'qy': np.array([0.2, 0.2, 0.0, 0.0, 0.0, 0.0])},
        'C4H2 + hv => C2H + C2H': {'qy': np.array([0.03, 0.03, 0.01, 0.01, 0.0, 0.0])},
        'C4H2 + hv => C2H2 + C2': {'qy': np.array([0.1, 0.1, 0.06, 0.06, 0.0, 0.0])},
        'C4H2 + hv => C4H2': {'qy': np.array([0.67, 0.67, 0.93, 0.93, 1.0, 1.0])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,vulcan), ('Wogan','VULCAN'),xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(wv2[0])], 'citations': ['Ferradaz2009']},
        {'nm-range': [float(wv2[0]), float(np.max(out['wv']))], 'citations': ['Smith1998']}
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