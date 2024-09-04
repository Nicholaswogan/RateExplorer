import utils
import numpy as np

def main():
    species = 'H2'
    out = utils.get_VULCAN(species)
    ld = utils.get_leiden(species)
    phid = utils.get_phidrates(species)

    out = utils.change_xs_to_other(out, ld)
    min_xs = np.min(out['wv'])
    out = utils.prepend_xs_of_other(out, phid)

    # Make plots
    utils.make_xs_plot(species, out, (ld, phid), ('Leiden','Phidrates'),xlim=(0,125))
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