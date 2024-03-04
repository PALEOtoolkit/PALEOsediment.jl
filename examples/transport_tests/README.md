# Sediment transport test cases

## Solid and solute transport

    include("PALEO_transport_sediment.jl")

Test solid and solute transport in sediment column with varying porosity.

## Multi-G particulate organic carbon transport and decay

    include("PALEO_transport_RCmultiG.jl")

Tests `ReactionRCmultiG` implementation of [Dale2015](@cite) discrete multi-G
representation of reactive-continuum particulate organic carbon. 

POC decay is added to DIC pool, includes d13C carbon isotopes.

Physical environment is from [Dale2015](@cite) 'shelf' and 'slope' cases.
