
import numpy as np
import yaml
import csv

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

def FormatReactions_main(data):
    
    order = ['reverse-reactions','atoms','species','particles','reactions','missing','reaction-notes']
    copy = data.copy()
    data.clear()
    for key in order:
        if key in copy.keys():
            data[key] = copy[key]
            
    # Atmos
    if 'atoms' in data:
        for i in range(len(data['atoms'])):
            data['atoms'][i] = flowmap(data['atoms'][i])
    
    # Species
    if 'species' in data:
        for i in range(len(data['species'])):
            
            if data['species'][i]['name'] == False:
                data['species'][i]['name'] = "NO"
            
            order = ['name', 'composition', 'condensate', 'thermo','note']
            copy = data['species'][i].copy()
            data['species'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['species'][i][key] = copy[key]
                        
            data['species'][i]['composition'] = flowmap(data['species'][i]['composition'])
            if 'thermo' in data['species'][i].keys():
                
                order = ['model', 'reference-pressure','temperature-ranges','data']
                copy = data['species'][i]['thermo'].copy()
                data['species'][i]['thermo'].clear()
                for key in order:
                    if key in copy.keys():
                        data['species'][i]['thermo'][key] = copy[key]
                    
                data['species'][i]['thermo']['temperature-ranges'] = blockseqtrue(data['species'][i]['thermo']['temperature-ranges'])
                
                data['species'][i]['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(data['species'][i]['thermo']['data'])]

    # Particles
    if 'particles' in data:
        for i in range(len(data['particles'])):
            data['particles'][i]['composition'] = flowmap(data['particles'][i]['composition'])
            if data['particles'][i]['formation'] == 'reaction':
                flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
                for key in flowstyle:
                    if key in data['particles'][i].keys():
                        data['particles'][i][key] = flowmap(data['particles'][i][key])
            elif data['particles'][i]['formation'] == 'saturation':
                flowstyle = ['parameters','vaporization','sublimation','super-critical']
                for key in flowstyle:
                    data['particles'][i]['saturation'][key] = flowmap(data['particles'][i]['saturation'][key])
            
    # Reactions
    if 'reactions' in data:
        for i in range(len(data['reactions'])):
            order = ['equation','type','rate-constant','rate-constants','low-P-rate-constant',
                     'high-P-rate-constant','duplicate','efficiencies','JPL','citation','location','checked', 'typo','note','ref']
            copy = data['reactions'][i].copy()
            data['reactions'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['reactions'][i][key] = copy[key]
                    
            flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies','citation','ref']
            for key in flowstyle:
                if key in data['reactions'][i].keys():
                    if isinstance(data['reactions'][i][key],dict):
                        data['reactions'][i][key] = flowmap(data['reactions'][i][key])

            if 'location' in data['reactions'][i]:
                if isinstance(data['reactions'][i]['location'],list):
                    data['reactions'][i]['location'] = blockseqtrue(data['reactions'][i]['location'])

            if 'rate-constants' in data['reactions'][i]:
                for j in range(len(data['reactions'][i]['rate-constants'])):
                    data['reactions'][i]['rate-constants'][j] = flowmap(data['reactions'][i]['rate-constants'][j])

    if 'reaction-notes' in data:
        for i in range(len(data['reaction-notes'])):
            order = ['reactions','note']
            copy = data['reaction-notes'][i].copy()
            data['reaction-notes'][i].clear()
            for key in order:
                if key in copy:
                    data['reaction-notes'][i][key] = copy[key]
            data['reaction-notes'][i]['reactions'] = blockseqtrue(data['reaction-notes'][i]['reactions'])
                    
    return data

def get_count(s):
    if s == 'O1D':
        s = 'O'
    elif s == 'N2D':
        s = 'N'
    
    # Use dictionary to store elements in insertion order
    mp = {}
 
    i = 0
    while i < len(s):
        count = 0
 
        # Convert the string element into character
        c = s[i]
 
        # Check string contains the Capital A - Z value.
        if c.isupper():
            a = " "
 
            # Convert the character element
            # into string element
            a += c
 
            j = i + 1
            while j < len(s):
                d = s[j]
 
                # Check string contains the small a-z value.
                if d.islower():
                    a += d
 
                    if a not in mp:
                        mp[a] = 1
                    else:
                        mp[a] += 1
                    count = 1
 
                # Check string contains the number value.
                elif d.isdigit():
                    k = int(d)
                    mp[a] = k
                    count = 1
                else:
                    i = j - 1
                    break
                j += 1
 
            if count == 0:
                if a not in mp:
                    mp[a] = 1
                else:
                    mp[a] += 1
 
        i += 1

    mp1 = {}
    for key in mp:
        mp1[key.strip()] = mp[key]
    return mp1

def compare2reactions(rx1, rx2):
    rx1 = rx1.replace('<=>','=>').replace('(','').replace(')','')
    rx2 = rx2.replace('<=>','=>').replace('(','').replace(')','')
    react1, prod1 = [sorted(a.replace(' ','').split('+')) for a in rx1.split('=>')]
    react1 = [a.lower() for a in react1]
    prod1 = [a.lower() for a in prod1]
    react2, prod2 = [sorted(a.replace(' ','').split('+')) for a in rx2.split('=>')]
    react2 = [a.lower() for a in react2]
    prod2 = [a.lower() for a in prod2]
    return all([react1 == react2, prod1 == prod2])

def compare2reactions_stoich(rx1, rx2):
    rx1 = rx1.replace('<=>','=>').replace('(','').replace(')','')
    rx2 = rx2.replace('<=>','=>').replace('(','').replace(')','')
    react1, prod1 = [sorted(a.replace(' ','').split('+')) for a in rx1.split('=>')]
    # react1 = [a.upper() for a in react1]
    # prod1 = [a.upper() for a in prod1]
    react2, prod2 = [sorted(a.replace(' ','').split('+')) for a in rx2.split('=>')]
    # react2 = [a.upper() for a in react2]
    # prod2 = [a.upper() for a in prod2]

    if len(react1) != len(react2):
        return False
    if len(prod1) != len(prod2):
        return False

    for i,r in enumerate(react1):
        if get_count(react1[i]) != get_count(react2[i]):
            return False
    for i,p in enumerate(prod1):
        if get_count(prod1[i]) != get_count(prod2[i]):
            return False
        
    return True

def compare2reactions_reactants(rx1, rx2):
    rx1 = rx1.replace('<=>','=>').replace('(','').replace(')','')
    rx2 = rx2.replace('<=>','=>').replace('(','').replace(')','')
    react1, prod1 = [sorted(a.replace(' ','').split('+')) for a in rx1.split('=>')]
    react1 = [a.lower() for a in react1]
    prod1 = [a.lower() for a in prod1]
    react2, prod2 = [sorted(a.replace(' ','').split('+')) for a in rx2.split('=>')]
    react2 = [a.lower() for a in react2]
    prod2 = [a.lower() for a in prod2]
    return all([react1 == react2])

def compare2rates(rate1,rate2):
    for key in rate1:
        if not np.isclose(rate1[key],rate2[key],rtol=1e-1,atol=1e-200):
            return False
    return True

def comparerates(rx1, rx2):
    keys = ['rate-constant','low-P-rate-constant','high-P-rate-constant']
    for key in keys:
        if key in rx1 and key not in rx2:
            return False
        if key not in rx1 and key in rx2:
            return False
    for key in keys:
        if key in rx1 and key in rx2:
            if not compare2rates(rx1[key],rx2[key]):
                return False
    return True

def comparerates_noA(rx1, rx2):
    keys = ['rate-constant','low-P-rate-constant','high-P-rate-constant']
    for key in keys:
        if key in rx1 and key not in rx2:
            return False
        if key not in rx1 and key in rx2:
            return False
    for key in keys:
        if key in rx1 and key in rx2 and key != 'A':
            if not compare2rates(rx1[key],rx2[key]):
                return False
    return True

def parse_latex_rate(rxn):
    
    rate_str = rxn.replace('$','').strip()
    rate = {'A': 0.0, 'b': 0.0, 'Ea': 0.0}
    if r'\!\times\!' in rate_str:
        tmp = rate_str.split(r'\!\times\!')
        a1 = tmp[0]
        ind = tmp[1].index('}')
        b1 = tmp[1][:ind+1].split('{')[1].split('}')[0]

        A = float(a1)*10**float(b1)

        b = 0.0
        if r'\left(' in rate_str:
            ind = rate_str.index(r'\right)^{')
            b = float(rate_str[ind:].split('}')[0].split('{')[1])

        Ea = 0.0
        if 'e^{' in rate_str:
            ind = rate_str.index('e^{')
            Ea = -float(rate_str[ind:].split('}')[0].split('{')[1].split('/')[0])

        A = A*(1/300)**b
        
        rate = {'A': float(A), 'b': float(b), 'Ea': float(Ea)}
    else:
        assert float(rate_str) == 0.0

    return rate

def latex_to_sp(a):
    return a.replace('_','').replace('$','').replace('(','').replace(')','').replace('^','')

def parse_latex_equation(reacts, prods):
    reacts1 = [latex_to_sp(a.strip()) for a in reacts.split('+')]
    prods1 = [latex_to_sp(a.strip()) for a in prods.split('+')]
    return ' + '.join(reacts1)+' <=> '+' + '.join(prods1)

def parse_latex_citation(a):
    res = a.replace(r'\\','').strip().split('%')[0].strip()
    if len(res) == 0:
        return 'Assumed'
    return res

def parse_reactions_latex():
    with open('all_reactions_half_auto2.tex','r',encoding="ISO-8859-1") as f:
        lines = f.readlines()
    lines = lines[32:861]

    rxns = []
    for i,line in enumerate(lines):
        if not line.strip().startswith(r'\refstepcounter{reaction}'):
            continue

        eq = {}
        
        cols = line.split('&')
        equation = parse_latex_equation(cols[1], cols[3])
        rate = parse_latex_rate(cols[4])
        citation = parse_latex_citation(cols[5])

        eq['equation'] = equation
        eq['rate-constant'] = rate
        eq['citation'] = citation

        if 'M' in equation:
            line1 = lines[i+1]
            if len(line1.strip()) > 0 and not line1.strip().startswith(r'\refstepcounter{reaction}'):
                cols = line1.split('&')
                rate1 = parse_latex_rate(cols[4])
                citation1 = parse_latex_citation(cols[5])
                del eq['rate-constant']
                del eq['citation']
                eq['type'] = 'falloff'
                eq['low-P-rate-constant'] = rate
                eq['high-P-rate-constant'] = rate1
                eq['citation'] = {'low-P': citation, 'high-P': citation1}
            else:
                eq['type'] = 'three-body'    
        rxns.append(eq)
    return rxns

def parse_reactions_latex2016():
    with open('Eridani_re_revised.tex','r',encoding="ISO-8859-1") as f:
        lines = f.readlines()
    lines = lines[1691:1790]

    k = 0
    rxns = []
    for i,line in enumerate(lines):
        if not line.strip().startswith(r'\refstepcounter{reaction}'):
            continue

        eq = {}
        
        cols = line.split('&')
        equation = parse_latex_equation(cols[1], cols[3])
        rate = parse_latex_rate(cols[4])
        citation = parse_latex_citation(cols[5])

        eq['equation'] = equation
        eq['rate-constant'] = rate
        eq['citation'] = citation

        if 'M' in equation:
            line1 = lines[i+1]
            if len(line1.strip()) > 0 and not line1.strip().startswith(r'\refstepcounter{reaction}'):
                cols = line1.split('&')
                rate1 = parse_latex_rate(cols[4])
                citation1 = parse_latex_citation(cols[5])
                del eq['rate-constant']
                del eq['citation']
                eq['type'] = 'falloff'
                eq['low-P-rate-constant'] = rate
                eq['high-P-rate-constant'] = rate1
                eq['citation'] = {'low-P': citation, 'high-P': citation1}
            else:
                eq['type'] = 'three-body'    
        rxns.append(eq)
        k += 1
    return rxns

def parse_citations_csv():
    with open('citations.csv', 'r') as file:
        csv_reader = csv.reader(file)
        rows = []
        for row in csv_reader:
            rows.append(row)
            # print(row)

    data = {
        'id': [a[0] for a in rows],
        'R1': [a[2] for a in rows],
        'R2': [a[3] for a in rows],
        'R3': [a[4] for a in rows],
        'P1': [a[5] for a in rows],
        'P2': [a[6] for a in rows],
        'P3': [a[7] for a in rows],
        'A': [a[9] for a in rows],
        'b': [a[10] for a in rows],
        'Ea': [a[11] for a in rows],
        'cite1': [a[14] for a in rows],
        'cite2': [a[15] for a in rows],
    }
    for i in range(len(data['A'])):
        try:
            data['A'][i] = float(data['A'][i])
        except ValueError:
            data['A'][i] = 0.0
        try:
            data['b'][i] = float(data['b'][i])
        except ValueError:
            data['b'][i] = 0.0
        try:
            data['Ea'][i] = float(data['Ea'][i])
        except ValueError:
            data['Ea'][i] = 0.0

    data1 = []
    for i in range(len(data['A'])):
        tmp = {key:data[key][i] for key in data}
        if tmp['R1'] == '':
            continue
        a = ' '.join([tmp['R1'],tmp['R2'],tmp['R3']])
        a = ' + '.join(a.split())
        b = ' '.join([tmp['P1'],tmp['P2'],tmp['P3']])
        b = ' + '.join(b.split())
        bb = b.split('+')
        if len(bb) == 1 and bb[0] != '' and bb[0].lower() != 'products':
            b += ' + M'
            a += ' + M'
        rxn = a + ' => '+b

        rate = {'A': tmp['A']*(1/300)**tmp['b'], 'b': tmp['b'], 'Ea': tmp['Ea']}

        rx = {'equation': rxn}
        rx['rate-constant'] = rate
        rx['citation'] = tmp['cite1']
        if tmp['cite2'] != '':
            rx['citation'] += '; '+tmp['cite2']
        if len(rx['citation']) == 0:
            rx['citation'] = 'Assumed CSV'
        rx['line'] = i + 1

        data1.append(rx)
 
    data = data1

    return data

def reactions_with_citations_from_csv(rxns, rxns_csv):
    rxns1 = []

    for rx in rxns:

        # Skip photolysis
        if 'hv' in rx['equation']:
            continue
        
        # Determine if falloff
        falloff = False
        if 'type' in rx:
            if rx['type'] == 'falloff':
                falloff = True

        found = False
        if not falloff:
            citation = None
            for rx1 in rxns_csv:
                if compare2reactions(rx['equation'],rx1['equation']) and compare2rates(rx['rate-constant'],rx1['rate-constant']):
                    citation = rx1['citation']
                    line = rx1['line']
                    found = True
                    break
        else:
            citation = {'low-P': None, 'high-P': None}
            for rx1 in rxns_csv:
                if compare2reactions(rx['equation'],rx1['equation']) and compare2rates(rx['low-P-rate-constant'],rx1['rate-constant']):
                    citation['low-P'] = rx1['citation']
                    line = rx1['line']
                    found = True
                    break

            for rx1 in rxns_csv:
                if compare2reactions(rx['equation'],rx1['equation']) and compare2rates(rx['high-P-rate-constant'],rx1['rate-constant']):
                    citation['high-P'] = rx1['citation']
                    line = rx1['line']
                    found = True
                    break

        # Skip reactions that are not found
        if not found:
            continue    
        
        rx['citation'] = citation
        rx['line'] = line
        rxns1.append(rx)

    return rxns1

def null_citation(falloff):
    if not falloff:
        return None
    else:
        return {'low-P': None, 'high-P': None}

def add_citations_to_zahnle():

    with open('zahnle_earth_v0.20.yaml','r') as f:
        sol = yaml.load(f,yaml.Loader)
    rxns = sol['reactions']

    rxns_latex = parse_reactions_latex()
    rxns_latex16 = parse_reactions_latex2016()
    rxns_csv1 = parse_citations_csv()
    rxns_csv = reactions_with_citations_from_csv(rxns, rxns_csv1)

    kk1 = 0
    kk2 = 0
    rxns_new = []
    for i,rx in enumerate(rxns):
        rx['checked'] = False
        
        falloff = False
        if 'type' in rx:
            if rx['type'] == 'photolysis':
                rxns_new.append(rx)
                continue
            elif rx['type'] == 'falloff':
                falloff = True

        # # Search latex
        # rx['location'] = []
        # found = False
        # for j,rx_l in enumerate(rxns_latex):
        #     if compare2reactions(rx['equation'],rx_l['equation']) and comparerates(rx, rx_l):
        #         rx['citation'] = rx_l['citation']
        #         rx['location'].append('latex R'+str(j+1))
        #         found = True
        #         break
        # if found:
        #     # If found, search for an equation match with latex doc
        #     for j,rx_c in enumerate(rxns_csv):
        #         if compare2reactions(rx['equation'],rx_c['equation']):
        #             rx['location'].append('equation-rate csv L'+str(rx_c['line']))
        #             break
        #     rxns_new.append(rx)
        #     continue

        # Search latex 2016
        # rx['location'] = []
        found = False
        for j,rx_l in enumerate(rxns_latex16):
            if compare2reactions(rx['equation'],rx_l['equation']) and comparerates(rx, rx_l):
                rx['citation'] = rx_l['citation']
                rx['location'] = 'latex16 R'+str(j+1)
                found = True
                break
        if found:
            rxns_new.append(rx)
            continue

        # return

        # # Search CSV
        # rx['location'] = []
        # found = False
        # for j,rx_c in enumerate(rxns_csv):
        #     if compare2reactions(rx['equation'],rx_c['equation']) and comparerates(rx, rx_c):
        #         rx['citation'] = rx_c['citation']
        #         rx['location'].append('csv L'+str(rx_c['line']))
        #         found = True
        #         break
        # if found:
        #     # If found, search for an equation match with latex doc
        #     for j,rx_l in enumerate(rxns_latex):
        #         if compare2reactions(rx['equation'],rx_l['equation']):
        #             rx['location'].append('equation-rate latex R'+str(j+1))
        #             break
        #     rxns_new.append(rx)
        #     continue

        kk1 += 1
        
        # rx['location'] = []
        # rx['citation'] = null_citation(falloff)
        # # Search latex for reactants + rate match
        # found = False
        # for j,rx_l in enumerate(rxns_latex):
        #     if compare2reactions_reactants(rx['equation'],rx_l['equation']) and comparerates_noA(rx, rx_l):
        #         rx['location'].append('reactants+rate latex R'+str(j+1))
        #         found = True
        #         break
        # if not found:
        #     # Search latex for equation match
        #     found = False
        #     for j,rx_l in enumerate(rxns_latex):
        #         if compare2reactions(rx['equation'],rx_l['equation']):
        #             rx['location'].append('equation-rate latex R'+str(j+1))
        #             found = True
        #             break
        #     if not found:
        #         # Search latex for reactants match
        #         found = False
        #         for j,rx_l in enumerate(rxns_latex):
        #             if compare2reactions_reactants(rx['equation'],rx_l['equation']):
        #                 rx['location'].append('reactants-rate latex R'+str(j+1))
        #                 found = True
        #                 break
        
        # # Search csv for reactants + rate match
        # found = False
        # for j,rx_c in enumerate(rxns_csv1):
        #     if compare2reactions_reactants(rx['equation'],rx_c['equation']) and comparerates_noA(rx, rx_c):
        #         rx['location'].append('reactants+rate csv L'+str(rx_c['line']))
        #         found = True
        #         break
        # if not found:
        #     # Search csv for equation match
        #     found = False
        #     for j,rx_c in enumerate(rxns_csv1):
        #         if compare2reactions(rx['equation'],rx_c['equation']):
        #             rx['location'].append('equation-rate csv L'+str(rx_c['line']))
        #             found = True
        #             break
        #     if not found:
        #         # Search csv for reactants match
        #         found = False
        #         for j,rx_c in enumerate(rxns_csv1):
        #             if compare2reactions_reactants(rx['equation'],rx_c['equation']):
        #                 rx['location'].append('reactants-rate csv L'+str(rx_c['line']))
        #                 found = True
        #                 break

        # if len(rx['location']) > 0:
        #     rxns_new.append(rx)
        #     continue

        kk2 += 1

        # No citation was found, so put in nulls
        rx['citation'] = null_citation(falloff)
        rx['location'] = 'None'
        rxns_new.append(rx)

    sol['reactions'] = rxns_new
    print(kk1, kk2)

    sol = FormatReactions_main(sol)
    with open('zahnle_earth_v0.20_cite1.yaml','w') as f:
        yaml.dump(sol, f, MyDumper,sort_keys=False,width=70)

if __name__ == '__main__':
    add_citations_to_zahnle()
    # parse_reactions_latex2016()

        








