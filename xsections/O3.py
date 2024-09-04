import utils
import numpy as np

def main():
    species = 'O3'
    out = utils.get_VULCAN(species)
    ld = utils.get_leiden(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    min_xs = np.min(out['wv'])

    # Use Leiden
    out = utils.change_xs_to_other(out,ld)

    # Prepend Phidrates, but assume it is all ionization
    other = phid
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds]*0,out['xsp'])

    # Adjust quantum yield
    out['ratios']['wv'] = np.append(out['ratios']['wv'],np.array([410.0, 411.0]))
    out['ratios']['O3 + hv => O + O2']['qy'] = np.append(out['ratios']['O3 + hv => O + O2']['qy'],np.array([out['ratios']['O3 + hv => O + O2']['qy'][-1], 1.0]))
    out['ratios']['O3 + hv => O1D + O2']['qy'] = np.append(out['ratios']['O3 + hv => O1D + O2']['qy'],np.array([out['ratios']['O3 + hv => O1D + O2']['qy'][-1], 0.0]))

    # Make plots
    utils.make_xs_plot(species, out, (phid,), ('Phidrates',),xlim=(0,1200))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(min_xs)], 'citations': ['Huebner2015']},
        {'nm-range': [float(min_xs), float(np.max(out['wv']))], 'citations': ['Heays2017']}
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), 220.0], 'citations': ['Huebner2015']},
        {'nm-range': [220.0, float(np.max(out['wv']))], 'citations': ['Matsumi2002']}
        ],
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()