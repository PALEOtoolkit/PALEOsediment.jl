# Sediment 1D reaction-transport FOAM evaluation

Configuration similar to FOAM test case from Zhao etal (2020) <https://dx.doi.org/10.1038/s41467-020-15673-3>


## Configuration with Corg, NO3, P, Mn, Fe, SO4, H2S, CH4 + minerals

    julia> include("PALEO_examples_sediment_FOAM.jl")

This configuration is similar to that of Zhao etal (2020) <https://dx.doi.org/10.1038/s41467-020-15673-3>,
with a 3m sediment column with 17 solute components, 18 solid components, 
and 12 Corg reactivity fractions approximating a reactive-continuum distribution.

The model output shows four columns, with a sensitivity test for Biotite input flux
showing the effect this has on pH.

This configuration should run to steady-state (at 1e5 yr) in ~1-2 minutes on a PC or laptop with 4 sediment columns,
with ~40 pseudo-transient-continuation iterations.

Differences from Zhao (2020):
- Organic matter uses fixed Corg:P stoichiometry for all reactivity fractions, 
  consequently P is high lower in the sediment.
- Phosphorus includes Fe co-precipitation, but only a minimal parameterisation of CFA
  (just a dependency on P concentration and F concentration)
- Uses TAlk instead of total protons as the component corresponding to H+ or charge balance.
  Chemically this is equivalent, but the diffusive transport need not be.
  Transport that includes Coulombic coupling would be better (at which point pH and speciation
  follow from charge balance). 
- Adsorbed Fe contribution to TAlk is implemented directly by including budget and transport contributions.

## Minimal config with SO4, H2S, CH4 only

    julia> include("PALEO_examples_sediment_FOAM.jl")

Minimal configuration with no N, Mn, Fe, pH and speciation, or minerals
