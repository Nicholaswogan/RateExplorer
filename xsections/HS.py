import utils
import numpy as np

def main():
    species = 'HS'
    out = utils.get_VULCAN(species)
    ld = utils.get_leiden('SH')
    # phid = utils.get_phidrates('SH')
    wogan = utils.get_wogan(species)

    out = utils.change_xs_to_other(out, ld)

    # Make sure absorption does not exceed photolysis and ionization
    inds = np.where(out['xsa'] < out['xsp'] + out['xsi'])
    out['xsa'][inds] = out['xsp'][inds] + out['xsi'][inds]

    # Make plots
    utils.make_xs_plot(species, out, (ld,wogan,), ('Leiden','Wogan',))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Heays2017']}
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