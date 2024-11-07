# -*- coding: utf-8 -*-
module SedimentReservoirs

# import Infiltrator # Julia debugger

import PALEOboxes as PB
using PALEOboxes.DocStrings


"""
    ReactionSedSolidReservoir, ReactionSedSolidReservoirTotal

A single (vector) reservoir (state variable) representing a solid phase biogeochemical tracer in a sediment

Parameter `molar_volume` defines density, used to calculate `R_volume_frac` (fraction of solid phase volume)
and add a contribution to `volume_change_sms (m^3 yr-1)` calculated from `R_sms`

Creates `R_conc` (mol m-3) and `R_conc_sms` (mol m-3 yr-1) as state variable and source-sink provided to the solver, 
and calculates `R` (mol).

Accumulation of fluxes into _sms is split into `R_sms` and `R_trspt_sms` to allow calculation of advective velocity from `volume_change_sms`:
- biogeochemical transformations and diffusive transport (bioturbation) should be added to `R_sms` (mol yr-1) 
(these are included in `volume_change_sms`)
- advective transport should be added to `R_trspt_sms` (not included in `volume_change_sms`).
- `R_conc_sms` adds up contributions from `R_sms` and `R_trspt_sms` for the solver.

In addition:
- if parameter `field_data <: AbstractIsotopeScalar` (eg `IsotopeLinear`), a Property `R_delta` is created.
- `ReactionSedSolidReservoirTotal` also calculates the Domain total `R_total` (units mol), eg to check budgets.

Local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Initialisation
Initial value is set using `variable_attributes:` in the configuration file, using `R:initial_value` and `R:initial_delta` 
(for 'state_conc=false') or `R_conc:initial_value` and `R_conc:initial_delta` (for a 'state_conc=true').
(NB: even if `initial_value` is set on `R`, it *sets concentration in `mol m-3`*.)

Transport is defined by attributes `:advect`, `:vertical_movement` (m d-1) set on the concentration variable `R_conc`. Optical
extinction is defined by the `:specific_light_extinction` (m^2 mol-1) attribute set on the concentration variable `R_conc`.

# Example configuration in .yaml file
                reservoir_Corg:  # sediment Corg
                    class: ReactionSedSolidReservoirTotal       # include _total (mol)
                    parameters:
                        molar_volume:                   12e-6 # m^3 / mol M  (Cx 30 g mol,  density 2.5 g / cm3  = 30 / 2.5 * 1e-6 )
                    variable_links:
                        R*:                             Corg*
                        volume:                         volume_solid
                    variable_attributes:
                        R_conc:vphase:                  VP_Solid                        
                        R:initial_value:                0.0  # concentration m-3 solid phase
                        R:norm_value:                   8.3e4  # mol m-3 solid  (~ 1 / molar_volume)

# See also
ReactionReservoir

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedSolidReservoir{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

        PB.ParBool("total", false,
            description="true to calculate R_total"),

        PB.ParDouble("molar_volume", 1e3, units="m^3 mol-1",
            description="molar volume"),

        PB.ParDouble("limit_delta_conc", 0.0, units="mol m-3",
            description="**EXPERIMENTAL** attempt to limit delta for low/-ve concentrations (0.0 to disable)"),
    )
end

abstract type ReactionSedSolidReservoirTotal <: PB.AbstractReaction end
function PB.create_reaction(::Type{ReactionSedSolidReservoirTotal}, base::PB.ReactionBase)
    rj = ReactionSedSolidReservoir(base=base)
    PB.setvalueanddefault!(rj.pars.total, true)
    PB.setfrozen!(rj.pars.total)
    return rj
end


function PB.register_methods!(rj::ReactionSedSolidReservoir)

    PB.setfrozen!(rj.pars.field_data)
    PB.setfrozen!(rj.pars.total)
    PB.setfrozen!(rj.pars.molar_volume)

    R_attributes=(
        :field_data=>rj.pars.field_data[],
        :calc_total=>rj.pars.total[],
    )
    R_conc_attributes = (
        :field_data=>rj.pars.field_data[], 
        :advect=>true,
        :vertical_movement=>0.0,
        :specific_light_extinction=>0.0,
        :vphase=>PB.VP_Undefined,
        :diffusivity_speciesname=>"",
        :diffusivity=>missing,
        :charge=>missing,
        :gamma=>missing,
    )

    volume  = PB.VarDep("volume",   "m3",       "cell volume (or cell phase volume eg for a sediment with solid and liquid phases)")


    R = PB.VarProp("R",        "mol",      "vector reservoir"; attributes=R_attributes)
    R_volume_frac = PB.VarProp("R_volume_frac",        "",      "fraction of solid phase volume")              
    R_conc = PB.VarStateExplicit("R_conc",   "mol m-3",  "concentration"; attributes=R_conc_attributes)
    PB.add_method_setup_initialvalue_vars_default!(rj, [R_conc])
    
    do_vars = PB.VariableReaction[R, R_volume_frac, R_conc, volume,]

    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(do_vars, PB.VarProp("R_delta", "per mil", "isotopic composition"))
        if rj.pars.limit_delta_conc[] > 0.0
            @warn "$(PB.fullname(rj)) using experimental limit_delta_conc $(rj.pars.limit_delta_conc[]) mol m-3"
        end
    end

    # calculate conc, volume etc
    PB.add_method_do!(rj, do_reactionsedsolidreservoir, (PB.VarList_namedtuple(do_vars),))

    # accumulate contribution to solid phase volume change
    R_sms = PB.VarTarget("R_sms",    "mol yr-1", "vector reservoir biogeochemical source-sinks",
        attributes=(:field_data=>rj.pars.field_data[], ),
    )
    volume_change_sms = PB.VarContrib("volume_change_sms",           "m^3 yr-1", "change in volume due to production and consumption of solid phase components")
    PB.add_method_do!(rj, do_reactionsedsolidreservoir_volume_change, (PB.VarList_namedtuple([R_sms, volume_change_sms]),))


    # calculate conc_sms
    R_conc_sms = PB.VarDeriv("R_conc_sms",    "mol m-3 yr-1", "vector reservoir source-sinks",
        attributes=(:field_data=>rj.pars.field_data[], ),
    )    
    R_trspt_sms = PB.VarTarget("R_trspt_sms",    "mol yr-1", "vector reservoir transport source-sinks",
        attributes=(:field_data=>rj.pars.field_data[], ),
    )
    PB.add_method_do!(rj, do_reactionsedsolidreservoirconc_sms, (PB.VarList_namedtuple([volume, R_conc_sms, PB.VarDep(R_sms), R_trspt_sms]),))
    
    if rj.pars.total[]
        PB.add_method_do_totals_default!(rj, [R])
    end
    PB.setfrozen!(rj.pars.total)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_reactionsedsolidreservoir(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        vars.R[i]  = vars.R_conc[i]*vars.volume[i]
        vars.R_volume_frac[i] = vars.R[i] * pars.molar_volume[] / vars.volume[i]
    end
    
    if hasfield(typeof(vars), :R_delta)
        limit_value = pars.limit_delta_conc[]
        if limit_value > 0.0
            # norm_value = PB.get_attribute(rj.var_R, :norm_value)::Float64
            # limit_value = 1e-6  # mol m-3

            @inbounds for i in cellrange.indices
                vars.R_delta[i]  = PB.get_delta_limit(vars.R[i], limit_value*vars.volume[i], 100.0)
            end
        else
            @inbounds for i in cellrange.indices
                vars.R_delta[i]  = PB.get_delta(vars.R[i])
            end
        end
    end

    return nothing
end

function do_reactionsedsolidreservoir_volume_change(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        vars.volume_change_sms[i]  += vars.R_sms[i]*pars.molar_volume[]
    end

    return nothing
end


function do_reactionsedsolidreservoirconc_sms(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        vars.R_conc_sms[i]  += (vars.R_sms[i] + vars.R_trspt_sms[i])/vars.volume[i]
    end

    return nothing
end



end # module
