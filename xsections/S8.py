import utils
import numpy as np

def main():
    species = 'S8'
    wv, xs = np.loadtxt('data/Bass1953_S8.txt').T
    out = {}
    out['wv'] = wv/10
    out['xsp'] = xs*0
    out['xsa'] = xs
    out['xsi'] = xs*0  

    wogan = utils.get_wogan(species)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,), ('Wogan',),xlim=None)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['']},
        ],
    "notes": "For S8 rings."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()