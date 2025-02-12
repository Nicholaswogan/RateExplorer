import yaml
from photochem.utils._format import yaml, Loader, MyDumper, flowmap, FormatReactions_main

def parse_zahnle_thermofile(filename):
    "Parses Kevin's thermodynamic files"

    with open(filename,'r') as f:
        lines = f.readlines()
    
    species = []
    for line in lines:
        species.append(line.split()[0])
    species = list(set(species))
    
    inds = []
    kinds = []
    for sp in species:
        ind = []
        kind = None
        for i,line in enumerate(lines):
            tmp = line.split()
            if tmp[0] == sp:
                ind.append(i)
                if len(tmp) == 12:
                    kind = 'normal'
                elif len(tmp) == 4:
                    kind = 'reference'
                    break
                else:
                    raise Exception(line)
        inds.append(ind)
        kinds.append(kind)
    
    # Do normals
    data = {}
    for i,sp in enumerate(species):
        if kinds[i] == 'normal':
            lbs = []
            ubs = []
            coeffs = []
            for ind in inds[i]:
                line = lines[ind]
                tmp = line.split()
                coeff = [float(a) for a in tmp[5:]]
                coeffs.append(coeff)
                lbs.append(1e3*float(tmp[1]))
                ubs.append(1e3*float(tmp[2]))
                enthalpy = float(tmp[3])
                entropy = float(tmp[4])
            entry = {}
            entry['temperature-ranges'] = lbs + [ubs[-1]]
            for i in range(len(entry['temperature-ranges'])-1):
                if entry['temperature-ranges'][i+1] <= entry['temperature-ranges'][i]:
                    raise Exception(sp)
            entry['data'] = coeffs
            entry['enthalpy'] = enthalpy
            entry['entropy'] = entropy
            entry['reference-species'] = None
            data[sp] = entry
    
    # Do reference
    for i,sp in enumerate(species):
        if kinds[i] == 'reference':
            line = lines[inds[i][0]]
            sp_name, enthalpy, entropy, sp_ref = line.split()
            enthalpy = float(enthalpy)
            entropy = float(entropy)
    
            thermo = data[sp_ref].copy()
            deltaH = enthalpy - thermo['enthalpy'] # enthalpy - reference enthalpy
            deltaS = entropy - thermo['entropy']
            
            coeffs = []
            for coeff in thermo['data']:
                coeff1 = coeff.copy()
                tmp = coeff1[5] + deltaH
                tmp = float('%.8e'%tmp)
                coeff1[5] = tmp
                tmp = coeff1[6] + deltaS
                tmp = float('%.8e'%tmp)
                coeff1[6] = tmp
                coeffs.append(coeff1)
            entry = {}
            entry['data'] = coeffs
            entry['enthalpy'] = enthalpy
            entry['entropy'] = entropy
            entry['temperature-ranges'] = thermo['temperature-ranges'].copy()
            entry['reference-species'] = sp_ref
            data[sp] = entry
    
    return data

def species_from_zahnlethermo(species_list, zahnle_thermofile):

    data = parse_zahnle_thermofile(zahnle_thermofile)
    
    with open('composition.yaml','r') as f:
        comp = yaml.load(f,Loader=yaml.Loader)
    
    out = []
    for sp in species_list:
        entry = {}
        entry['name'] = sp
        entry['composition'] = comp[sp]
        entry['thermo'] = {}
        entry['thermo']['model'] = 'Shomate'
        entry['thermo']['temperature-ranges'] = data[sp]['temperature-ranges']
        entry['thermo']['data'] = data[sp]['data']
        if data[sp]['reference-species'] is None:
            entry['note'] = 'From the NIST database'
        else:
            tmp = data[sp]['reference-species']
            entry['note'] = 'Estimated from thermodynamic data at 298 K and species '+tmp
        out.append(entry)
    
    species = {}
    species['species'] = out
    return species

def main():
    with open('composition.yaml','r') as f:
        comp = yaml.load(f,Loader=Loader)

    zahnle_thermofile = 'thermodata121_wHe_Tlim.rx'
    species_list = [key for key in comp]
    species = species_from_zahnlethermo(species_list, zahnle_thermofile)
    species = FormatReactions_main(species)

    with open('thermodata121.yaml','w') as f:
        yaml.dump(species,f,sort_keys=False,Dumper=MyDumper,width=70)

if __name__ == '__main__':
    main()
