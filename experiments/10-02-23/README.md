
# 10-02-23

This repository contains updates for the reactions H + H2CO (+ M) <=> CH3O (+ M) and H + H2CO (+ M) <=> H2COH (+ M), along with the thermodynamics of CH3O and H2COH. The new data is all based on [Xu et al. (2015)](https://doi.org/10.1021/acs.jpca.5b00553). The new thermodynamics is reflected in Kevin Zahnle's `thermodata121.rx`.

I make comparisons to Shami's (i.e., Shang-Min Tsai) rates and thermodynamic data. The rates come from the `NCHO_photo_network.yaml` in VULCAN: https://github.com/exoclime/VULCAN/blob/acddf3b6490bc587c90876dfe406a5755ecbeff8/thermo/NCHO_photo_network.txt

In summary, here are the new thermodynamics and rates that should be used

**New thermodynamics**

```yaml
- name: CH3O
  composition: {C: 1, H: 3, O: 1}
  thermo:
    model: Shomate
    temperature-ranges: [0.0, 1200.0, 6000.0]
    data:
    - [5.193767, 93.23249, -44.85457, 7.882279, 0.551175, 25.5409, 217.5163]
    - [71.35268, 6.174497, -1.19109, 0.079564, -15.58917, -25.7327, 277.368]
  note: Estimated from thermodynamic data at 298 K and species H2CO

- name: H2COH
  composition: {H: 3, C: 1, O: 1}
  thermo:
    model: Shomate
    temperature-ranges: [0.0, 3000.0, 6000.0]
    data:
    - [20.559, 95.615, -31.306, 3.5141, 0.025109, -19.0, 241.9]
    - [111.5, 3.0, 0.0, 0.0, 0.0, -85.7, 310.6]
  note: Estimated from thermodynamic data at 298 K and species CH3OH
```

**New rates**

```yaml
# Fit to Xu et al. (2015)
- equation: H + H2CO (+ M) <=> CH3O (+ M)
  type: falloff
  low-P-rate-constant: {A: 1.13e-20, b: -3.973, Ea: 3281.8}
  high-P-rate-constant: {A: 8.96e-12, b: -0.205, Ea: 1877.0}

# Fit to Xu et al. (2015)
- equation: H + H2CO (+ M) <=> H2COH (+ M)
  type: falloff
  low-P-rate-constant: {A: 4.98e-28, b: -1.625, Ea: 3032.0}
  high-P-rate-constant: {A: 5.80e-09, b: -1.095, Ea: 3074.9}
```

