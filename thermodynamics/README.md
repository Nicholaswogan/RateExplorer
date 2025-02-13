# Zahnle data

This folder holds Kevin Zahnle's thermodynamic data and reactions. Details are also documented in this repository: https://github.com/Nicholaswogan/ImpactAtmosphere/blob/main/ImpactAtmosphere/data/notes.md

# 10/2/23

I have added thermodata121.rx on the basis of the following email of Kevin Zahnle on 9/29/23. We chose new enthalpies of formation for CH3O and H2COH based on [Xu et al. (2015)](https://doi.org/10.1021/acs.jpca.5b00553). Xu gives enthalpies at 0 K, but we can convert them to 298 K following this website: https://cccbdb.nist.gov/enthform2.asp

# 2/12/25

I added thermodata121_wHe_Tlim.rx which adds He, and changes the T limits of the polynomial fits to be consistent with what I'm currently using in photochem v0.6.3.

I also created thermodata122.yaml, building off of thermodata121_wHe_Tlim.rx. For all species except He and H2, I created a new set of polynomials for 10 K < T < 298 K. The polynomials are designed to smoothly vary at the 298 K transition, and also at low T tend toward a heat capacity that only has translational and rotational degrees of freedom:

$$c_p = \left(\frac{5}{2}\right)R + \left(\frac{N_\mathrm{rot}}{2}\right)R$$

Here, $N_\mathrm{rot}$, is the number of rotational degrees of freedom. The file `rotational_dof.yaml` has the number of rotational degrees of freedom for each species.

Now, the thermodynamics of all species are reasonable down to hopefully ~10 K.