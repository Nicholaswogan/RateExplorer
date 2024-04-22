# 4-18-24

Here, I put together a new rate for the `HO2 + HO2 <=> H2O2 + O2` reaction based on [Klippenstein et al. (2022)](https://doi.org/10.1016/j.combustflame.2021.111975) and [Kircher and Sander (1984)](https://doi.org/10.1021/j150654a029). Klippenstein is an Ab-initio study that gives the rate at 400 to 2000 K, while Kircher and Sander estimate the rate below 400 K with experimental data. Note that the reaction appears to have pressure dependence at low temperatures (< 400 K), according to Kircher and Sander. Here, I ignore this pressure dependence because it is relatively minor between 0 and 1 bar.

The script `fit_rate.py` constructs a new rate using Kircher and Sander for T < 400 K and Klippenstein for T > 400 K. The script also produces a plot that compares the new rate to Kircher and Sander and Klippenstein.

The new rate below is valid for 150 K to 2000 K. The rate is likely too slow for low temperatures (T < 300 K) and high pressures (P > 10 bar), but these conditions are relatively unlikely in many planetary atmospheres.

```yaml
- equation: HO2 + HO2 <=> H2O2 + O2
  rate-constant: {A: 3.19e-19, b: 2.01, Ea: -1274.3}
```