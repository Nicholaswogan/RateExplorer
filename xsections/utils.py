import h5py
import numpy as np
from urllib.request import urlopen
from io import BytesIO
from zipfile import ZipFile
import os
import yaml
from matplotlib import pyplot as plt

# Yaml stuff
class flowmap( dict ): pass
def flowmap_rep(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=True)
class blockseqtrue( list ): pass
def blockseqtrue_rep(dumper, data):
    return dumper.represent_sequence( u'tag:yaml.org,2002:seq', data, flow_style=True )
yaml.add_representer(blockseqtrue, blockseqtrue_rep)
yaml.add_representer(flowmap, flowmap_rep)
class MyDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()

def download_and_unzip(url, extract_to='.'):
    http_response = urlopen(url)
    zipfile = ZipFile(BytesIO(http_response.read()))
    zipfile.extractall(path=extract_to)

# Lieden dataset
if not os.path.isdir('cross_sections'):
    download_and_unzip('https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/all_cross_sections_text_continuum_0.1nm.zip')

# PhotoData, including phidrates
if not os.path.isdir('PhotoData'):
    download_and_unzip('https://github.com/Nicholaswogan/PhotoData/archive/5c259fae902072081ba2d90390d5cb6e02134264.zip')
    os.rename('PhotoData-5c259fae902072081ba2d90390d5cb6e02134264','PhotoData')
from PhotoData.PhotoData import phidrates, MPI_Mainz

# VULCAN
if not os.path.isdir('VULCAN'):
    download_and_unzip('https://github.com/exoclime/VULCAN/archive/f3d7291d69b356a38f18d70a39c41e143eb85cee.zip')
    os.rename('VULCAN-f3d7291d69b356a38f18d70a39c41e143eb85cee','VULCAN')

# Zahnle data
if not os.path.isdir('photochem_clima_data'):
    download_and_unzip('https://github.com/nicholaswogan/photochem_clima_data/archive/a1f8cd1324c438d3543775855eabb059f5ba1c5b.zip')
    os.rename('photochem_clima_data-a1f8cd1324c438d3543775855eabb059f5ba1c5b','photochem_clima_data')
with open('photochem_clima_data/reaction_mechanisms/zahnle_earth.yaml','r') as f:
    ZAHNLE_DATA = yaml.load(f,yaml.Loader)

def get_leiden(species):
    filename = 'cross_sections/'+species+'/'+species+'_0.1nm.txt'
    with open(filename,'r') as f:
        lines = f.readlines()
    line = lines[2]
    labels = line.split()[1:]
    if labels == ['wavelength', 'photoabsorption', 'photodissociation', 'photoionisation']:
        wv,xsa,xsp,xsi = np.loadtxt(filename,skiprows=3).T
    elif labels == ['wavelength', 'photoabsorption', 'photoionisation']:
        wv,xsa,xsi = np.loadtxt(filename,skiprows=3).T
        xsp = np.zeros(wv.shape[0])
    elif labels == ['wavelength', 'photoabsorption', 'photodissociation']:
        wv,xsa,xsp = np.loadtxt(filename,skiprows=3).T
        xsi = np.zeros(wv.shape[0])
    else:
        raise Exception('Problem reading '+filename)
    out = {}
    out['wv'] = wv
    out['xsa'] = xsa
    out['xsp'] =  xsp
    out['xsi'] = xsi
    return out

def get_mpimainz(species, best='single longest'):
    mpi = MPI_Mainz(species)
    mpi.get_data()
    mpi.find_best_data(best=best)
    out = {}
    out['wv'] = mpi.best_data['wavelength']
    out['xsa'] = mpi.best_data['cross section']
    out['xsp'] = mpi.best_data['cross section']
    out['xsi'] = mpi.best_data['cross section']*0.0

    return out

def get_phidrates(species):

    phi = phidrates(species)
    xsp = np.zeros(phi.neutral.data['wavelength'].shape[0])
    for i,a in enumerate(phi.neutral.branches):
        if a != 'SO2band': # the one exception
            xsp += phi.neutral.data[a]
    xsp *= phi.neutral.data['cross section']

    xsi = np.zeros(phi.neutral.data['wavelength'].shape[0])
    for i,a in enumerate(phi.ion.branches):
        xsi += phi.ion.data[a]
    xsi *= phi.ion.data['cross section']

    xsa = phi.neutral.data['cross section']

    out = {}
    out['wv'] = phi.neutral.data['wavelength']
    out['xsa'] = xsa
    out['xsp'] =  xsp
    out['xsi'] = xsi

    return out

def get_phidrates_yields(species):
    phi = phidrates(species)
    xsp = np.zeros(phi.neutral.data['wavelength'].shape)
    for branch in phi.neutral.branches:
        xsp += phi.neutral.data[branch]
    xsp = phi.neutral.data['cross section']*xsp
    inds = np.where(xsp > 0)
    
    qys = {'wv':phi.neutral.data['wavelength'][inds]}
    for branch in phi.neutral.branches:
        qys[branch] = (phi.neutral.data[branch][inds]*phi.neutral.data['cross section'][inds])/xsp[inds]

    return qys

def get_zahnle_data(species):

    with open('Zahnle_Kevin_data/SW_neutrals_for_stays.DAT','r') as f:
        lines = f.readlines()
    line = lines[0]
    labels = line.split()[4:]
    out = {}
    for key in labels:
        out[key] = np.empty((len(lines)-2)*2)
    out['wv'] = np.empty((len(lines)-2)*2)
    
    for i,line in enumerate(lines[2:]):
        ii = i*2
        tmp = [float(a) for a in line.split()]
        out['wv'][ii] = tmp[1]
        out['wv'][ii+1] = tmp[2]
        for j,key in enumerate(labels):
            out[key][ii] = tmp[j+4]
            out[key][ii+1] = tmp[j+4]

    with open('Zahnle_Kevin_data/SW_ions_for_stays.DAT','r') as f:
        lines = f.readlines()
    line = lines[0]
    labels = line.split()[4:]
    out1 = {}
    for key in labels:
        out1[key] = np.empty((len(lines)-2)*2)
    out1['wv'] = np.empty((len(lines)-2)*2)
    
    for i,line in enumerate(lines[2:]):
        ii = i*2
        tmp = [float(a) for a in line.split()]
        out1['wv'][ii] = tmp[1]
        out1['wv'][ii+1] = tmp[2]
        for j,key in enumerate(labels):
            out1[key][ii] = tmp[j+4]
            out1[key][ii+1] = tmp[j+4]
            
    assert np.all(out['wv'] == out1['wv'])
    
    out2 = {}
    out2['wv'] = out['wv']/10
    out2['wv'][0] = 0.1
    
    out2['xsp'] = np.zeros(len(out['wv']))
    if species in out:
        out2['xsp'] = out[species]

    out2['xsi'] = np.zeros(len(out['wv']))
    if species in out1:
        out2['xsi'] = out1[species]

    out2['xsa'] = out2['xsp'] + out2['xsi']

    if species not in out and species not in out1:
        raise Exception()

    return out2

def get_wogan(species):
    wv, xsp = np.loadtxt('photochem_clima_data/xsections/'+species+'/'+species+'_xs.txt',skiprows=2).T
    xsa = xsp.copy()
    xsi = np.zeros(xsa.shape[0])

    out = {}
    out['wv'] = wv
    out['xsa'] = xsa
    out['xsp'] =  xsp
    out['xsi'] = xsi

    ratio_data = {}
    qyfiles = [a for a in os.listdir('photochem_clima_data/xsections/'+species+'/') if 'hv' in a]
    for i,qy in enumerate(qyfiles):
        rx = qy.strip('.txt').replace('=','=>').replace('_',' ')
        wv, qy = np.loadtxt('photochem_clima_data/xsections/'+species+'/'+qy,skiprows=2).T

        if i == 0:
            ratio_data['wv'] = wv.copy()
        else:
            assert np.all(ratio_data['wv'] == wv)

        ratio_data[rx] = {}
        ratio_data[rx]['qy'] = qy
        ratio_data[rx]['new'] = False

    out['ratios'] = ratio_data
    out['missing'] = []

    return out

def minmaxarray(arr):
    return np.array([np.min(arr),np.max(arr)])

def compare2reactions(rx1, rx2):
    rx1 = rx1.replace('<=>','=>').replace('(','').replace(')','')
    rx2 = rx2.replace('<=>','=>').replace('(','').replace(')','')
    react1, prod1 = [sorted(a.replace(' ','').split('+')) for a in rx1.split('=>')]
    react2, prod2 = [sorted(a.replace(' ','').split('+')) for a in rx2.split('=>')]
    return all([react1 == react2, prod1 == prod2])

def rename_as_zahnle(rx, data):
    zahnle_photos = [a['equation'] for a in data['reactions'] if 'hv' in a['equation']]
    for rx1 in zahnle_photos:
        if compare2reactions(rx1, rx):
            return True,rx1
    return False,rx

def rename_all_as_zahnle(ratios):
    ratios_new = {'wv': ratios['wv']}
    for key in ratios:
        if key != 'wv':
            found, rx = rename_as_zahnle(key, ZAHNLE_DATA)
            new = not found
            ratios_new[rx] = ratios[key]
            ratios_new[rx]['new'] = new
    return ratios_new            

def get_zahnle_branches(sp, data):
    zahnle_photos = [a['equation'] for a in data['reactions'] if 'hv' in a['equation']]
    zbranches = []
    for rx in zahnle_photos:
        sp1 = rx.split('=>')[0].split('+')[0].strip()
        if sp1 == sp:
            zbranches.append(rx)
    return zbranches
            
def get_VULCAN(species):
    
    folder = 'VULCAN/thermo/'
    
    # Vulcan names to Photochem  names
    v2p = {
        'O_1':'O1D',
        'N_2D':'N2D',
        'CH2_1':'1CH2',
        'SH':'HS',
        'COS':'OCS'
    }
    vspecies = species
    for key in v2p:
        if v2p[key] == species:
            vspecies = key
            break
    
    # Download
    wv, xsa, xsp, xsi = np.loadtxt(folder+'photo_cross/'+vspecies+'/'+vspecies+'_cross.csv',skiprows=1,delimiter=',').T
    
    # Branching ratios
    ratios = np.genfromtxt(folder+'photo_cross/'+vspecies+'/'+vspecies+'_branch.csv',dtype=float,delimiter=',',skip_header=1, names = True)
    
    branches = {}
    with open(folder+'SNCHO_full_photo_network.txt') as f:
        lines = f.readlines()
    
    for i,line in enumerate(lines):
        if line.startswith('# photo disscoiation'):
            ind = i
            break
    for i,line in enumerate(lines[ind+2:-1]):
        react1 = line.split(']')[1].split()[0].strip()
        br = int(line.split(']')[1].split()[1].strip())
        rx = line.split('[')[1].split(']')[0]
        tmp = [a.strip() for a in rx.split('->')]
        prod = [a.strip() for a in tmp[1].split('+')]
        prod1 = []
        for p in prod:
            if p in v2p:
                prod1.append(v2p[p])
            else:
                prod1.append(p)
        prod = prod1
        react = tmp[0]
        assert react == react1
        if react in v2p:
            react = v2p[react]
    
        rx = react+' + hv => '+' + '.join(prod)
    
        if species == react:
            branches[br] = rename_as_zahnle(rx, ZAHNLE_DATA)

    ratio_data = {}
    ratio_data['wv'] = ratios['lambda']
    for i in range(len(ratios.dtype.names)-1):
        j = i + 1
        name = 'br_ratio_'+str(j)
        ratio_data[branches[j][1]] = {}
        ratio_data[branches[j][1]]['qy'] = ratios[name]
        ratio_data[branches[j][1]]['new'] = not branches[j][0]

    zmissing = get_missing_zahnle(species, ratio_data)
    
    out = {
        'wv':wv,
        'xsa':xsa,
        'xsp':xsp,
        'xsi':xsi,
        'ratios':ratio_data,
        'missing':zmissing
    }
    return out

def get_missing_zahnle(species, ratio_data):
    zbranches = get_zahnle_branches(species, ZAHNLE_DATA)
    zmissing = []
    for zbranch in zbranches:
        if zbranch not in ratio_data:
            zmissing.append(zbranch)
    return zmissing

def change_xs_to_other(out, other):
    out['wv'] = other['wv']
    out['xsa'] = other['xsa']
    out['xsp'] =  other['xsp']
    out['xsi'] = other['xsi']
    return out 

def prepend_xs_of_other(out, other, photo_must_be_zero=True):
    inds = np.where(other['wv']<out['wv'][0])
    out['wv'] = np.append(other['wv'][inds],out['wv'])
    out['xsa'] = np.append(other['xsa'][inds],out['xsa'])
    out['xsi'] = np.append(other['xsi'][inds],out['xsi'])
    out['xsp'] = np.append(other['xsp'][inds],out['xsp'])
    if photo_must_be_zero:
        assert np.all(other['xsp'][inds] == 0) # make sure photolysis is zero
    return out

def check_xs_and_qy(out):

    # Check for tiny numbers
    minval = 1e100
    for key in ['xsa','xsp','xsi']:
        if len(out[key][out[key] > 0]) > 0:
            minval = np.minimum(np.min(out[key][out[key] > 0]),minval)
    assert minval > 1e-45
    
    # Check to make sure that QYs sum to 1
    if 'ratios' in out:
        qyt = np.zeros(out['ratios']['wv'].shape[0])
        for ratio in out['ratios']:
            if ratio != 'wv': 
                qyt += out['ratios'][ratio]['qy']
        assert np.all(np.isclose(qyt ,1,rtol=1e-3))

    # Make sure absorption is never less than the sum of ionoization and photolysis
    # To within a tolerance
    assert np.all(out['xsa']*(1+1e-3) >= out['xsp'] + out['xsi'])

    # Must be sorted
    assert np.all(out['wv'][:-1] <= out['wv'][1:])

def make_h5_from_dict(species, out):

    check_xs_and_qy(out)
    
    outdir='results/'
    with h5py.File(outdir+species+'.h5','w') as f:
        dset = f.create_dataset("wavelengths", out['wv'].shape, 'f')
        dset[:] = out['wv']

        dset = f.create_dataset("photoabsorption", out['xsa'].shape, 'f')
        dset[:] = out['xsa']

        dset = f.create_dataset("photodissociation", out['xsp'].shape, 'f')
        dset[:] = out['xsp']

        dset = f.create_dataset("photoionisation", out['xsi'].shape, 'f')
        dset[:] = out['xsi']

        if 'ratios' in out:
            grp = f.create_group("photodissociation-qy")

            dset = grp.create_dataset("wavelengths", out['ratios']['wv'].shape, 'f')
            dset[:] = out['ratios']['wv']
            
            for ratio in out['ratios']:
                if ratio != 'wv':
                    dset = grp.create_dataset(ratio, out['ratios'][ratio]['qy'].shape, 'f')
                    dset[:] = out['ratios'][ratio]['qy']   

def make_xs_plot(species, out, compare=(), compare_labels=(), xlim=None, ylim=(1e-30,1e-15), save=True):
    plt.rcParams.update({'font.size': 14})
    fig,axs = plt.subplots(1,3,figsize=[15,4],sharey=True)

    keys = ['xsp','xsi','xsa']
    for i,ax in enumerate(axs):

        for j,tmp in enumerate(compare):
            ax.plot(tmp['wv'], tmp[keys[i]], alpha=0.6,lw=2,ls='-', label=compare_labels[j])

        ax.plot(out['wv'], out[keys[i]], alpha=0.6,c='k',lw=2,ls='--', label='XS')

    labels = ['Photolysis','Ionization','Tot. absorption']
    for i,ax in enumerate(axs):
        ax.set_yscale('log')
        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.grid(alpha=0.4)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_title(labels[i])

    ax = axs[0]
    ax.set_ylabel('Cross section (cm$^2$/molecule)')
    ax = axs[1]
    ax.legend(ncol=3,bbox_to_anchor=(0.5, 1.1), loc='lower center')
    
    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    if save:
        plt.savefig('figures/'+species+'_xs.pdf',bbox_inches='tight')
        plt.close()

def make_qy_plot(species, out, save=True, print_branches=True):
    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(1,1,figsize=[6,4],sharey=True)

    qyt = np.zeros(out['ratios']['wv'].shape[0])
    for ratio in out['ratios']:
        if ratio != 'wv': 
            ax.plot(out['ratios']['wv'],out['ratios'][ratio]['qy'],marker='o',label=ratio)
            qyt += out['ratios'][ratio]['qy']

    # Make sure QY sum to 1
    assert np.all(np.isclose(qyt ,1,rtol=1e-3))

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Quantum Yield')
    ax.legend(ncol=1,bbox_to_anchor=(0.5, 1.1), loc='lower center')

    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    if save:
        plt.savefig('figures/'+species+'_qy.pdf',bbox_inches='tight')
        plt.close()

    if print_branches:

        news = []
        for ratio in out['ratios']:
            if ratio != 'wv':
                if out['ratios'][ratio]['new']:
                    news.append(ratio)

        if len(news) == 0:
            print('New Branches: None')
        else:
            print('New Branches:')
            for new in news:
                print('- '+new)

        missings = []
        for branch in out['missing']:
            missings.append(branch)

        if len(missings) == 0:
            print('Branches missing from Zahnle: None')
        else:
            print('Branches missing from Zahnle:')
            for missing in missings:
                print('- '+missing)

def format_citation(species, citation):
    for i in range(len(citation[species]['xsections'])):
        citation[species]['xsections'][i]['nm-range'] = blockseqtrue(citation[species]['xsections'][i]['nm-range'])
        citation[species]['xsections'][i]['citations'] = blockseqtrue(citation[species]['xsections'][i]['citations'])
        citation[species]['xsections'][i] = flowmap(citation[species]['xsections'][i])
    if 'photodissociation-qy' in citation[species]:
        for i in range(len(citation[species]['photodissociation-qy'])):
            citation[species]['photodissociation-qy'][i]['nm-range'] = blockseqtrue(citation[species]['photodissociation-qy'][i]['nm-range'])
            citation[species]['photodissociation-qy'][i]['citations'] = blockseqtrue(citation[species]['photodissociation-qy'][i]['citations'])
            citation[species]['photodissociation-qy'][i] = flowmap(citation[species]['photodissociation-qy'][i])
    return citation

def save_citation(species, citation):
    citation = format_citation(species, citation)
    with open('results/'+species+'.yaml','w') as f:
        yaml.dump(citation, f, Dumper=yaml.Dumper, sort_keys=False)

