import utils
import numpy as np

def main():
    species = 'SO2'
    out = utils.get_VULCAN(species)
    ld = utils.get_leiden(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    out = utils.change_xs_to_other(out, ld)

    # Make plots
    utils.make_xs_plot(species, out, (ld, phid, wogan), ('Leiden','Phidrates', 'Wogan'),xlim=(0,500))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
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