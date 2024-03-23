# Sediment 1D reaction-transport configuration with C, O, N, P, S, Fe, Mn

## Dale (2015) configuration

    julia> include("PALEO_examples_sediment_NNFeMn.jl")

Configuration is similar to shelf and slope test cases in 
Dale (2015) GBC <https://dx.doi.org/10.1002/2014GB005017>, with:
- a discrete (multi-G) approximation to POC reactive-continuum distribution
- NO2 and NO3 oxidants
- Mn oxides with two reactivity fractions
- Fe oxides with three reactivity fractions

There are four sediment columns: oxic and anoxic 'shelf', oxic and anoxic 'slope'

## Minimal model for P diagenesis

    julia> include("PALEO_examples_sediment_NNFeMnP.jl")

Includes a minimal model for P diagenesis including:
  - P coprecipitation on FeHR, FeMR iron oxides
  - CFA formation

See:
- Slomp (1996) <https://dx.doi.org/10.1357/0022240963213745>
- Reed (2011) <https://dx.doi.org/10.4319/lo.2011.56.3.1075>
- Dale (2016) <https://dx.doi.org/10.1016/j.gca.2016.05.046>



    
