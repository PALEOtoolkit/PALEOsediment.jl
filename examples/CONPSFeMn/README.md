# Sediment 1D reaction-transport configuration with C, O, N, P, S, Fe, Mn

    julia> include("PALEO_examples_sediment_NNFeMn.jl")

Configuration is similar to shelf and slope test cases in 
Dale (2015) GBC <https://dx.doi.org/10.1002/2014GB005017>, with:
- a discrete (multi-G) approximation to POC reactive-continuum distribution
- NO2 and NO3 oxidants
- Mn oxides with two reactivity fractions
- Fe oxides with three reactivity fractions

There are four sediment columns: oxic and anoxic 'shelf', oxic and anoxic 'slope'

Differences from Dale (2015):

- Fe-S system:
    - solves speciation between FeII, H2S and FeSaq
      (Rickard & Luther 2007 <https://dx.doi.org/10.1021/cr0503658>)
      with FeSm preciptiation based on FeSaq solubility threshold.
    - Pyrite formation only includes 'Berzelius' pathway
      (FeSm + H2S -> FeS2pyr + H2)
    - No S0 

- Includes a minimal model for P diagenesis including:
    - P coprecipitation on FeHR, FeMR iron oxides
    - CFA formation

    See:
    - Slomp (1996) <https://dx.doi.org/10.1357/0022240963213745>
    - Reed (2011) <https://dx.doi.org/10.4319/lo.2011.56.3.1075>
    - Dale (2016) <https://dx.doi.org/10.1016/j.gca.2016.05.046>
    - Zhao (2020) <https://dx.doi.org/10.1038/s41467-020-15673-3>


    
