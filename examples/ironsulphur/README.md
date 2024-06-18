# Sediment iron-sulphur chemistry examples

## Implementation
Fe-S chemistry is implemented by three PALEO reactions ReactionFeSaq, ReactionFeSm, ReactionPyrH2S,
see PALEOaqchem documentation for details

#### ReactionFeSaq: equilibrium chemistry of FeII, H2S, FeS system

Represents two equilibrium reactions as algebraic constraints to
define primary species [Fe2++] and [HS-] in terms of:
- total [S-II] = [HS-] + [H2S] + [FeSaq]
- total [FeII] = [Fe++] + [FeSaq] 
given a provided fixed value of [H+]

    H2S        <--> H+ + HS-         eqb. const. K_H2S
    FeSaq + H+ <--> Fe++ + HS-       eqb. const K2_FeSaq

#### ReactionFeSm: FeS precipitation/dissolution

Represents 

    [FeSaq] <--> [FeSm]

as a fast precipitation dissolution reaction

#### ReactionPyrH2S: Pyrite formation by the "H2S" (Berzelius) mechanism

    FeSm + H2S -> FeS2pyr + H2

at rate 

    R_Pyr_H2S * [FeSm] * [H2S]

# Sediment configurations and experiments

NB: these are illustrative model test configurations only, and are not scientifically validated !!

Physical environment and organic matter input from PALEOexamples\src\sediment\boudreau1996\README.md,
from [Boudreau1996](@cite) shelf and shelf/slope cases.

Fe reactions with FeHR, FeMR, FePR iron oxide phases from Dale (2015), Van Cappellen & Wang (1996), with Fe-S system with
explicit treatment of FeS from  Rickard and Luther (2007), Lenton & Daines (2017), van de Velde (2021).

TODO 
- not checked rate constants carefully
- physical configuration and organic matter input is illustrative only
- no Mn

## Oxygen, sulphur, iron, with iron-sulphur interaction and pyrite, iron oxide burial

    julia> include("PALEO_examples_sediment_Fe_pyr_x15.jl")

Set of 15 columns, all with same shelf environment, Corg input with two reactivity fractions (a combination of the high-reactivity fraction from the shelf/slope case, and the low reactivity fraction from the rise case), excess Fe input:
- Columns 1-5 with [SO4] 28 mM (high modern value), decreasing values of [O2], with bioturbation
- Columns 6-10 with [SO4] 28 mM (high modern value), decreasing values of [O2], no bioturbation
- Columns 11-15 with zero [O2], no bioturbation decreasing values of [SO4]

Parameters for Fe cycle modified to create a high Fe test case, and increase FeS2pyr oxidation
- FeHR input ~0.6 mol m-2 yr-1 (comparable to Van Cappellen & Wang 1996 test case, ~x10 higher than a standard modern shelf value?)
- decrease FeMR limiting concentration for remineralization by x2 vs Van Cappellen & Wang 1996 value (NB: the "Monod" inhibition here is not the same as the linear form they used)
- increase FeS2pyr oxidation rate constant x10

Demonstrates (cf van de Velde & Meysman (2016)):
- Pyrite oxidation with modern [O2] and bioturbated sediment resulting in (2*FeS2)/(S reduction) of ~0.1, increases once [O2] drops below ~100 uM,
  cf Canfield & Farquar (2009)
- (non bioturbated oxic sediment with modern [O2] stays mostly oxic with little S reduction hence low pyrite burial)
- Zero oxygen (and non bioturbated) sediment shows high SR, (2*FeS2)/(S reduction) of ~0.5 hence high pyrite burial until  [SO4] drops to ~100uM
  (down a factor of two at ~1 Mm), presumably set by the remin [SO4] dependence (inhibition below 1 mM), cf Habicht (2002)

## Test case with pyrite burial

    julia> include("PALEO_examples_sediment_Fe_pyr.jl")

Boudreau (1996) test cases with Corg, O2, SO4/H2S, P, Fe-S and pyrite burial.
Three columns: shelf/slope,  rise, rise (no bioturbation)

POC stoichiometry Corg, P only (no N) to simplify test of alkalinity budget, ie
check alkalinity flux from SO4 and FeII solute fluxes = soluteflux_TAlk (as solid phase inputs and
output make no contribution to alkalinity budget).

    julia> PALEOmodel.get_array(paleorun.output, "fluxOceanfloor.soluteflux_SmIIaqtot", (cell=1, tmodel=1e12)).values
    0.00021099944013347013

    julia> PALEOmodel.get_array(paleorun.output, "fluxOceanfloor.soluteflux_SO4", (cell=1, tmodel=1e12)).values
    -0.0745449641168272  # net sulphate input balanced by FeS / FeS2 burial

    julia> PALEOmodel.get_array(paleorun.output, "fluxOceanfloor.soluteflux_FeIIaqtot", (cell=1, tmodel=1e12)).values
    0.0012053409002682764

    julia> PALEOmodel.get_array(paleorun.output, "fluxOceanfloor.soluteflux_TAlk", (cell=1, tmodel=1e12)).values
    0.15150061203704088

    julia> -2*-0.0745449641168272 + 2*0.0012053409002682764
    0.15150061003419096 # check alkalinity flux from SO4 and FeII fluxes = soluteflux_TAlk

## Test case with Fe only, no Fe-S or pyrite burial

    julia> include("PALEO_examples_sediment_Fe.jl")

Boudreau (1996) test cases with Corg, O2, SO4/H2S, P, Fe, no Fe-S or pyrite
Three columns: shelf/slope,  rise, rise (no bioturbation)

# References

- Canfield & Farquar (2009) PNAS <https://dx.doi.org/10.1073/pnas.0902037106>
- Habicht etal (2002) Science <https://dx.doi.org/10.1126/science.1078265>
- Poulton & Canfield (2011) Elements <https://dx.doi.org/10.2113/gselements.7.2.107>
- Rickard (2006) GCA <https://dx.doi.org/10.1016/j.gca.2006.02.029>
- Rickard and Luther (2007) Chemical Reviews <https://dx.doi.org/10.1021/cr0503658>
- Lenton and Daines (2017) Ann. Rev. Mar. Sci. <https://dx.doi.org/10.1146/annurev-marine-010816-060521>
- van de Velde & Meysman (2016) Aquatic Geochem. <https://dx.doi.org/10.1007/s10498-016-9301-7>
- van de Velde etal (2021) GMD <https://dx.doi.org/10.5194/gmd-14-2713-2021>
- Dale etal (2015) GBC <https://dx.doi.org/10.1002/2014GB005017>
- Van Cappellen & Wang Am J Sci (1996) <https://dx.doi.org/10.2475/ajs.296.3.197>
