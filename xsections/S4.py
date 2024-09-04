import utils
import numpy as np

def get_xs():
    wv, xsp = np.loadtxt('data/Billmers1991_S4.txt').T
    out = {}
    out['wv'] = wv
    out['xsp'] = xsp
    out['xsa'] = xsp
    out['xsi'] = xsp*0
    return out

def main():
    species = 'S4'
    out = get_xs()

    vulcan = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    ratios = {
        'wv': np.array([np.min(out['wv']),np.max(out['wv'])]),
        'S4 + hv => S3 + S': {'qy': np.array([0.5, 0.5]), 'new': False},
        'S4 + hv => S2 + S2': {'qy': np.array([0.5, 0.5]), 'new': True}
    }
    missing = []
    out['ratios'] = ratios
    out['missing'] = missing

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, wogan,), ('VULCAN','Wogan',))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])),float(np.max(out['wv']))], 'citations': ['KellerRudek2013','Billmers1991']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])),float(np.max(out['wv']))], 'citations': ['Assumed']}
        ],
    'notes': "Here, I assume even split between S3 + S and S2 + S2 chanel."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()