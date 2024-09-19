import utils
import numpy as np

def main():
    species = 'N'
    out = utils.get_leiden(species)
    phi = utils.get_phidrates(species)

    min_xs = np.min(out['wv'])
    out = utils.prepend_xs_of_other(out, phi)

    # Make plots
    utils.make_xs_plot(species, out, (phi,), ('Phidrates',), xlim=(0,300))

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()