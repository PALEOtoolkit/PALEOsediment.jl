# Equations, units, and PALEO variable name conventions

Conservation equations for a species `C` in PALEO are defined in terms of moles of `C` per cell, where a cell `volume (m^3)` is partitioned into solute volume `volume_solute = volume * phi (m^3)` and a solid volume `volume_solid = volume * phi_solid (m^3)` where `phi` is the porosity and `phi_solid = 1 - phi`.

The conservation equation for a species `C` is then:

    dC/dt = C_sms   (mol yr-1)

A solute concentration is defined as moles per m^3 of solute:

    C_conc = C / volume_solute  (mol m-3)

and a solid concentration as moles per m^3 of solid:

    C_conc = C /  volume_solid (mol m-3)

## Reservoir configuration

For the usual case where `C` is defined by a `ReactionReservoir` or `ReactionReservoirTotal`, the yaml file configuration for a solute reservoir should include:

                    variable_links:
                        volume:                         volume_solute   # relink volume so C_conc is solute concentration
                    variable_attributes: 
                        R_conc:vphase:                  VP_Solute      # label phase for transport etc

and the yaml file for a solid reservoir should include:

                    variable_links:
                        volume:                         volume_solid   # relink volume so C_conc is solid-phase concentration
                    variable_attributes: 
                        R_conc:vphase:                  VP_Solid      # label phase for transport etc

## Rate laws for reactions

PALEO units for rate constants for both solute and solid-phase first-order rate constants are `(mol m-3)-1 yr-1`.  This follows the convention in [VanCappellen1996a](@cite) (see Table 3 and equations 72 and 76).

#### Reactions between solute species

For a hypothetical reaction between solutes A, B, C

    A + B -> C

with first-order kinetics hence reaction rate

    R_abc = k_abc * A_conc * B_conc   (mol (m^3 solute)-1 yr-1)

`k_abc` has units `(mol m-3)-1 yr-1`, where concentration units are `mol (m^3 solute)-1`.  `A_sms` and similar contributions to tracer conservation equations are then:

    A_sms = - R_abc * volume_solute  (mol yr-1)

and similar for `B_sms`, `C_sms`.

This means that PALEO Reactions implementing eg secondary redox between solute species usually will require the yaml configuration to include:

                    variable_links:
                        volume:                         volume_solute   # relink volume so volume refers to solute volume

#### Reactions between a solid phase species and solute phase species

For a hypothetical reaction between solid-phase species `D` and solute `B`

    D + B -> C

with first-order kinetics hence reaction rate

    R_dbc = k_dbc * D_conc * B_conc   (mol (m^3 solid)-1 yr-1)

`k_dbc` (still) has units `(mol m-3)-1 yr-1`, where concentration units for `D` are `mol (m^3 solid)-1` and for B are `mol (m^3 solute)-1`.  `D_sms` and similar contributions to tracer conservation equations are then:

    D_sms = - R_dbc * volume_solid  (mol yr-1)

and similar for `B_sms`, `C_sms`.

This means that PALEO Reactions implementing eg secondary redox between a solid and solute species usually will require the yaml configuration to include:

                    variable_links:
                        volume:                         volume_solid   # relink volume so volume refers to solid phase volume

## Unit conversions for rate laws

    [mol/l]-1 [yr-1]    = [mol m-3] [yr-1]
    1 M-1 yr-1          = 1e-3  (mol m-3)-1 yr-1
    1 M-0.5 yr-1        = 3.16e-2 (mol m-3)^-0.5 yr-1

    1 M-1 h-1           = 8.766  (mol m-3)-1 yr-1
    1 M-0.5 h-1         = 277.2  (mol m-3)^-0.5 yr-1
