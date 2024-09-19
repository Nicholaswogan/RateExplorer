import utils
import numpy as np

def main():
    species = 'C'
    out = utils.get_leiden(species)
    phi = utils.get_phidrates('C3P')

    # Make plots
    utils.make_xs_plot(species, out, (phi,), ('Phidrates',), xlim=(0,300))

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']},
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()