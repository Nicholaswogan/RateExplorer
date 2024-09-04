import utils
import numpy as np

def main():
    species = 'CH3CHO'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)
    wogan = utils.get_wogan(species)

    qys = utils.get_phidrates_yields('CH3CHO')
    qys_new = {}
    qys_new['wv'] = qys['wv']
    qys_new['CH3CHO + hv => CH4 + CO'] = {'qy':qys['CH4/CO']}
    qys_new['CH3CHO + hv => CH3 + HCO'] = {'qy': qys['CH3/HCO'] + qys['tCH3/CHO']}
    ratios = qys_new
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (vulcan, phid, wogan), ('VULCAN','Phidrates', 'Wogan'),xlim=(0,400))
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