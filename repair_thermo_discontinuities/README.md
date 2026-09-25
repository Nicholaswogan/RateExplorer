# Repairing Zahnle Earth thermodynamic discontinuities

## 1. Why this repair is needed

`zahnle_earth.yaml` uses adjacent Shomate fits to describe each gas species over
different temperature ranges. At a boundary, those fits should give the same heat
capacity, enthalpy, entropy, and Gibbs energy. Many did not. The largest original
Gibbs-energy jump was $+106.434\,\mathrm{kJ\,mol^{-1}}$ for C3H6 at
$1400\,\mathrm{K}$; nine related species had jumps near
$-10.038\,\mathrm{kJ\,mol^{-1}}$ at $1500\,\mathrm{K}$.

The large jumps were mainly in enthalpy: ten joins exceeded
$1\,\mathrm{kJ\,mol^{-1}}$, including one near $100\,\mathrm{kJ\,mol^{-1}}$ and
nine near $10\,\mathrm{kJ\,mol^{-1}}$. Heat capacity and entropy also had
discontinuities, but they were much smaller. The largest heat-capacity jump was
$0.325\,\mathrm{J\,mol^{-1}\,K^{-1}}$ (less than 1% of the value at that join),
and only one entropy jump exceeded $1\,\mathrm{J\,mol^{-1}\,K^{-1}}$.
The large Gibbs-energy jumps were therefore driven mainly by enthalpy.

Photochem uses Gibbs energy to calculate temperature-dependent equilibrium factors and
reverse reaction rates. A jump therefore makes those rates change abruptly as an
atmospheric layer crosses a fit boundary. That disrupted integration in the bundled
TOI-1231 b case. This repair removes the jumps; it does not claim to replace the
underlying thermodynamic measurements.

## 2. Reproduce the repair

From the RateExplorer repository root:

```sh
cd repair_thermo_discontinuities
conda env create -f environment.yaml
conda run -n repair python run_repair.py
conda run -n repair python photochem_case/run_case.py
```

The environment pins Photochem **0.9.0**, `photochem_clima_data` **0.3.2**, and the core
Python libraries by version without platform-specific builds. `run_repair.py` retrieves
the data package at tag `v0.3.2` (commit `66789aba287a4eb8a846e3913ce346c0ceb5215f`) and
reads its original files from Git. It writes the two final files here:

- `zahnle_earth.yaml`: repaired gas thermodynamics and the unchanged reactions and saturation models.
- `condensate_thermo.yaml`: regenerated condensate thermodynamics.

The fetched source, numerical checks, and four comparison figures are under the ignored
`results/` directory. The repair does not modify the fetched checkout. The
`photochem_case/` folder holds its code, inputs, generated case YAMLs, and run log; its
generated files are also ignored. NASA9 comparisons fetch SHA-256-checked [NASA CEA
data](https://github.com/nasa/cea/blob/3f4441d28a02fccbb140e1a028d9902390981389/data/thermo.inp)
and the [Burcat/Goos/Ruscic database](https://respecth.elte.hu/burcat/NEWNASA.TXT).

## 3. Gas thermodynamics repair

For every join, `match_shomate_joins.py` keeps the earlier Shomate segment fixed and
changes only three coefficients in the next segment: `A` to match heat capacity, `F` to
match enthalpy, and `G` to match entropy. Since $G=H-TS$, it then
matches too. Coefficients `B` through `E`, temperature ranges, species, particles, and
reactions are unchanged. The correction proceeds through all later segments of each
species.

Enthalpy had the largest original discontinuities. C3H6 had a
$100\,\mathrm{kJ\,mol^{-1}}$ enthalpy jump at $1400\,\mathrm{K}$; its repair
produced the largest change along any gas curve, $127.1\,\mathrm{kJ\,mol^{-1}}$
in Gibbs energy. Nine other species had enthalpy jumps of about
$10\,\mathrm{kJ\,mol^{-1}}$ at $1500\,\mathrm{K}$: CH2CO, HCNOH, C2H2OH,
CH3CO, CH3O2, NH2CO, CH2N2, CH2CN, and CH3CN. All other enthalpy jumps were
below $1\,\mathrm{kJ\,mol^{-1}}$; the heat-capacity and entropy jumps were also
small. These major jumps are already present in the older
[`thermodata120_2-10-2021.rx`](../thermodynamics/thermodata120_2-10-2021.rx)
data. C3H6 and CH2CO have the same discontinuities in that file, and the other
eight species inherit CH2CO's fits with fixed enthalpy and
entropy offsets. As far as this repository can establish, our YAML conversion
and low-temperature extension did not introduce them. The repair makes all 181
joins continuous across the 96 species with multiple segments. The largest
remaining Gibbs-energy jump is
$2.33\times10^{-9}\,\mathrm{J\,mol^{-1}}$.

NASA9 supports the enthalpy repairs, which address the main problem: the
high-temperature enthalpy RMS disagreement improves for all six species with
unambiguous NASA9 matches. Gibbs-energy RMS worsens for five of them because
the original enthalpy error partly canceled a pre-existing entropy difference
from NASA9. The entropy changes very little in this repair, and its contribution
to $G=H-TS$ grows with temperature.

## 4. Condensate thermodynamics repair

The 13 particle saturation models in `zahnle_earth.yaml` stay unchanged: they were
fitted to vapor-pressure data and do not depend on the gas repair. The condensate fits
*do* depend on gas Gibbs energy, so `derive_condensates.py` regenerates them from the
repaired gas curves and those same saturation models. It fits condensate `A`, `F`, and
`G` over the original solid, liquid, and supercritical temperature windows, then adjusts
`F` at the triple and critical points to keep Gibbs energy continuous. Real
phase-transition jumps in enthalpy, entropy, and heat capacity remain.

This follows the same fitting objective and temperature windows as
[`fit_thermo`](../saturation_thermo/fitting.py): fit `A`, `F`, and `G` to log saturation
pressure, then align the branches at the triple and critical points. The
implementation differs. `fit_thermo` reads the older `thermodata121.yaml` gas fits and
separate `*_sat.yaml` files, uses a SciPy root solver, rounds to seven significant
figures, and writes one file per species. This script reads the repaired gas fits and
their unchanged, embedded saturation models, solves the linear least-squares fit
directly, adjusts `F` analytically, retains 15-digit coefficients, and writes the
combined `condensate_thermo.yaml`. Because the repaired gas Gibbs energy is continuous,
matching condensate Gibbs energy also matches saturation pressure at each boundary.

The largest remaining condensate Gibbs jump is
$4.41\times10^{-9}\,\mathrm{J\,mol^{-1}}$. Below the critical temperature, the
largest saturation-pressure error is $0.01993$ in $\log_{10}(P/\mathrm{bar})$
(about $4.7\%$), and the largest latent-heat error is
$2.35\,\mathrm{kJ\,mol^{-1}}$. These are fit residuals, not boundary jumps. The
supercritical branch is an extrapolation: its pressure error can reach $3.72$ in
$\log_{10}(P/\mathrm{bar})$ by $6000\,\mathrm{K}$, so it should not be read as
validated high-temperature phase equilibrium.

## 5. Photochem 0.9.0 validation

`photochem_case/run_case.py` regenerates the repaired files and case-specific
mechanisms, then runs the TOI-1231 b diagnostic with Photochem 0.9.0. The case contains
its own climate and spectrum inputs. The diagnostic **converged after 1,348 robust
steps** at model time $1.787\times10^{9}\,\mathrm{s}$, with **zero solver errors**. Its detailed
record is `photochem_case/photochem_case_run.jsonl`.

This shows that the repaired files work in the previously troublesome integration. It
does not by itself establish that every high-temperature species fit or supercritical
condensate extrapolation is physically correct.
