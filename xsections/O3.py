import utils
import numpy as np

def get_ratios():
    T = 250 # we assume 250 K
    wv = np.linspace(306,328,50)
    lam = wv

    X1, X2, X3 = 304.225, 314.957, 310.737
    omega1, omega2, omega3 = 5.576, 6.601, 2.187
    A1, A2, A3 = 0.8036, 8.9061, 0.1192
    nu1, nu2 = 0, 825.518
    c = 0.0765

    R = 0.695
    q1 = np.exp(-nu1/(R*T))
    q2 = np.exp(-nu2/(R*T))
    qy = (q1/(q1 + q2))*A1*np.exp(-((X1 - lam)/omega1)**4) + (q2/(q1 + q2))*A2*(T/300)**2*np.exp(-((X2 - lam)/omega2)**2) \
    + A3*(T/300)**1.5*np.exp(-((X3 - lam)/omega3)**2) + c

    # Recommends 0.9 from 220 to 305, and we take Phidrates below 220
    wv = np.append([219,220,305],wv)
    qy = np.append([1,0.9,0.9],qy)

    # Recommends 0.9 from 220 to 305
    wv = np.append(wv,[410.0,411.0])
    qy = np.append(qy,[qy[-1],0])

    ratios = {
        'wv': wv,
        'O3 + hv => O + O2': {'qy': 1-qy},
        'O3 + hv => O1D + O2': {'qy': qy},
    }

    return ratios

def main():
    species = 'O3'
    out = utils.get_leiden(species)

    vulcan = utils.get_VULCAN(species)
    phid = utils.get_phidrates(species)

    # Prepend Phidrates, but assume it is all ionization
    min_xs = np.min(out['wv'])
    other = phid
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsa'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds]*0,out['xsp'])

    # Quantum yields
    ratios = get_ratios()
    ratios = utils.rename_all_as_zahnle(ratios)
    out['ratios'] = ratios
    out['missing'] = utils.get_missing_zahnle(species, ratios)

    # Make plots
    utils.make_xs_plot(species, out, (phid,vulcan), ('Phidrates','VULCAN'),xlim=(0,1200))
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
    'notes': "I assume 250 K yields from Matsumi et al (2002) between 306 and 328 nm."
    }
    citation = {species: tmp}
    utils.save_citation(species, citation)

    # Save
    utils.make_h5_from_dict(species, out)

if __name__ == '__main__':
    main()