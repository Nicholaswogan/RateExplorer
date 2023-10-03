import utils

class XuCH3ORate:
    'H + H2CO (+ M) <=> CH3O (+ M)'
    def __init__(self):
        self.rate = {}
        self.rate['1atm_lowT'] = {'A': 2.32e-10, 'b': -1.22, 'Ea': 1813.2}
        self.rate['1atm_highT'] = {'A': 3.10e8, 'b': -6.79, 'Ea': 5573.9}
        self.rate['100atm_lowT'] = {'A': 7.12e-14, 'b': 0.45, 'Ea': 1587.1}
        self.rate['100atm_highT'] = {'A': 3.8e21, 'b': -10.15, 'Ea': 9000.0}

    def rate_1atm(self, T):
        if T < 800:
            tmp = self.rate['1atm_lowT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        else:
            tmp = self.rate['1atm_highT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k

    def rate_100atm(self, T):
        if T < 800:
            tmp = self.rate['100atm_lowT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        else:
            tmp = self.rate['100atm_highT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k
    
class XuCH3OReverseRate:
    'CH3O (+ M) <=> H + H2CO (+ M)'
    def __init__(self):
        self.rate = {}
        self.rate['1atm'] = {'A': 3.17e24, 'b': -4.25, 'Ea': 13104.9}
        self.rate['100atm'] = {'A': 8.5e28, 'b': -5.09, 'Ea': 14384.5}

    def rate_1atm(self, T):
        tmp = self.rate['1atm']
        k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k

    def rate_100atm(self, T):
        tmp = self.rate['100atm']
        k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k

class XuH2COHRate:
    'H + H2CO (+ M) <=> H2COH (+ M)'
    def __init__(self):
        self.rate = {}
        self.rate['1atm_lowT'] = {'A': 7e-10, 'b': -1.4, 'Ea': 2612.5}
        self.rate['1atm_highT'] = {'A': 3.14e7, 'b': -6.23, 'Ea': 7720.3}
        self.rate['100atm_lowT'] = {'A': 7.43e-21, 'b': 2.84, 'Ea': 3003.7}
        self.rate['100atm_highT'] = {'A': 6.09e15, 'b': -8.04, 'Ea': 10826.2}

    def rate_1atm(self, T):
        if T < 1000:
            tmp = self.rate['1atm_lowT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        else:
            tmp = self.rate['1atm_highT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k

    def rate_100atm(self, T):
        if T < 1000:
            tmp = self.rate['100atm_lowT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        else:
            tmp = self.rate['100atm_highT']
            k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k
    
class XuH2COHReverseRate:
    'H2COH (+ M) <=> H + H2CO (+ M)'
    def __init__(self):
        self.rate = {}
        self.rate['1atm'] = {'A': 4.52e34, 'b': -7.11, 'Ea': 22176.3}
        self.rate['100atm'] = {'A': 3.39e32, 'b': -5.88, 'Ea': 23371.8}

    def rate_1atm(self, T):
        tmp = self.rate['1atm']
        k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k

    def rate_100atm(self, T):
        tmp = self.rate['100atm']
        k = utils.arrhenius_rate(tmp['A'],tmp['b'],tmp['Ea'], T)
        return k