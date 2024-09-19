import utils
import numpy as np
from scipy.interpolate import interp1d

def main():
    species = 'H2SO4'

    wogan = utils.get_wogan(species)

    wv, xs = np.loadtxt('data/Lane2008_H2SO4.txt').T

    # Log linear extrapolation to 230 nm
    wv_new = np.linspace(wv[-1]*1.01,230,10)
    xs_new = 10.0**interp1d(wv,np.log10(xs),fill_value='extrapolate',bounds_error=False)(wv_new)
    wv = np.append(wv,wv_new)
    xs = np.append(xs,xs_new)

    out = {}
    out['wv'] = wv
    out['xsa'] = xs
    out['xsp'] = xs
    out['xsi'] = xs*0.0

    # Quantum yields
    ratios = {
        'wv': np.array([118.13, 181.29]),
        'H2SO4 + hv => SO3 + H2O': {'qy': np.array([1.0, 1.0])},
    }
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (wogan,), ('Wogan',), xlim=(0,400))
    utils.make_qy_plot(species, out)

    # Save citation
    tmp = {
    "xsections": [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Lane2008']},
        ],
    'photodissociation-qy': [
        {'nm-range': [float(np.min(out['wv'])), float(np.max(out['wv']))], 'citations': ['Zhang2012']}
        ],
    'notes': 'I extrapolated the Lane et al. (2008) cross section beyond 181 nm.'
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()