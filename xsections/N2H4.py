import utils
import numpy as np

def main():
    species = 'N2H4'
    out = utils.get_VULCAN(species)
    wogan = utils.get_wogan(species)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,), ('Wogan',),xlim=None, ylim=None)
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