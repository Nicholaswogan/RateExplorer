
# S + O2 <=> SO + O

I did a refit to the [Lu et al. (2004)](http://dx.doi.org/10.1063/1.1792611) rate for S + O2 <=> SO + O to make the extrapolation to low temperatures more sensible. This rate produces something reasonable from 50 K or 100 K to 3000 K

```yaml
- equation: S + O2 <=> SO + O
  rate-constant: {A: 3.21e-16, b: 1.35, Ea: -285.8}
```

# SO3 + H2O + H2O <=> H2SO4 + H2O

I checked, and it looks like the Lovejoy rate I already have matches Figure 5 in [Jayne et al. (1997)](https://doi.org/10.1021/jp972549z). Extrapolation to high or low T looks sensible. Also, the reverse of the Lovejoy rate doesn't look crazy for a wide range of T.

# CL + O3 <=> ClO + O2

I checked to make sure the rate extrapolated well to high and low T. Also, the rate matches more or less what is in JPL 19-5. I'm happy with the Atkinson rate.

# C2H4 + N2D <=> CH3CN + H

This reaction is from Lavvas et al. (2008). This reaction doesn't actually happen in reality, but I'm going to leave it anyway. The reactions temperature dependence isn't insane, and the rate is reasonable. Also, I need this reaction to reproduce CH3CN in Titan's atmosphere. Making CH3CN the correct way in Titan's atmosphere would take work because it requires adding molecules to the mechanism.
