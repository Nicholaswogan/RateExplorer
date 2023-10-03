
# 10-02-23

This repository contains updates for the following reactions
- H + H2CO (+ M) <=> CH3O (+ M)
- H + H2CO (+ M) <=> H2COH (+ M)
- CH3 + O <=> H2CO + H
- H2CO + H <=> HCO + H2

It also updates the thermodynamics for CH3O and H2COH. The new data is all based on [Xu et al. (2015)](https://doi.org/10.1021/acs.jpca.5b00553). The new thermodynamics is reflected in Kevin Zahnle's `thermodata121.rx`.

In summary, here are the new thermodynamics and rates that should be used

**New thermodynamics**

```yaml
- name: CH3O
  composition: {C: 1, H: 3, O: 1}
  thermo:
    model: Shomate
    temperature-ranges: [0.0, 1200.0, 6000.0]
    data:
    - [5.193767, 93.23249, -44.85457, 7.882279, 0.551175, 18.1409, 217.5163]
    - [71.35268, 6.174497, -1.19109, 0.079564, -15.58917, -33.1327, 277.368]
  note: Estimated from thermodynamic data at 298 K and species H2CO
- name: H2COH
  composition: {H: 3, C: 1, O: 1}
  thermo:
    model: Shomate
    temperature-ranges: [0.0, 3000.0, 6000.0]
    data:
    - [20.559, 95.615, -31.306, 3.5141, 0.025109, -25.3, 241.9]
    - [111.5, 3.0, 0.0, 0.0, 0.0, -92.0, 310.6]
  note: Estimated from thermodynamic data at 298 K and species CH3OH
```

**New rates**

```yaml
# Fit to Xu et al. (2015)
- equation: H + H2CO (+ M) <=> CH3O (+ M)
  type: falloff
  low-P-rate-constant: {A: 1.22e-23, b: -3.0, Ea: 2900}
  high-P-rate-constant: {A: 6.56e+03, b: -5, Ea: 4000}

# Fit to Xu et al. (2015)
- equation: H + H2CO (+ M) <=> H2COH (+ M)
  type: falloff
  low-P-rate-constant: {A: 2.82e-29, b: -1.2, Ea: 2900}
  high-P-rate-constant: {A: 3.0e-12, b: 0.0, Ea: 3500}

# Given by Xu et al. (2015)
- equation: CH3 + O <=> H2CO + H
  rate-constant: {A: 9.0e-11, b: 0.0, Ea: 0.0}
- equation: CH3 + O <=> HCO + H2
  rate-constant: {A: 6.0e-11, b: 0.0, Ea: 0.0}

# Given by Xu et al. (2015)
- equation: H2CO + H <=> HCO + H2
  rate-constant: {A: 2.28e-19, b: 2.65, Ea: 766.5}
```
