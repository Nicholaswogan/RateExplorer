import utils
import numpy as np

def main():
    species = 'Cl'
    leiden = utils.get_leiden(species)
    out = utils.get_phidrates(species)

    # Make plots
    utils.make_xs_plot(species, out, (leiden,), ('Leiden',), xlim=(0,300))

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Huebner2015']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()