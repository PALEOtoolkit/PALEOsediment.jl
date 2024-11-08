module SedimentPhys

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOsediment

import Printf

import Infiltrator # julia debugger

"""
    ReactionSedimentPhys

Calculate sediment physical parameters: porosity, advection velociy, salinity, temperature.

The sediment grid is defined by eg [`PALEOsediment.Sediment.SedimentGrid.ReactionSedimentGridn1D`](@ref).

Porosity is set to a prescribed time-independent functional form defined by `f_porosity` parameter.

Advection velocity may calculated in two ways:
- as a time-independent function of depth calculated from `oceanfloor.w_accum`, and assuming no volume changes in solid phase.
- calculated from solid-phase volume changes accumulated into `volume_change_sms` variable.

## Physical environment

Per-column environment is defined by `oceanfloor` variables.

Parameters for prescribed time-independent porosity are set by `oceanfloor.phi` and optionally `oceanfloor.phimin`

Prescribed accumulation rate is set by `oceanfloor.w_accum` (omit or set to zero to calculate 
accumulation rate from solid-phase volume changes).

Temperature and salinity are set by `oceanfloor.sal`, `oceanfloor.temp`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedimentPhys{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParString("f_porosity", "Const", allowed_values=["Const", "ExpAtten"],
            description="functional form for porosity vs depth"),
        PB.ParDouble("zpor",  0.1, units="m",
            description="lengthscale for porosity if f_porosity=ExpAtten"),
        PB.ParBool("w_solute", false,
            description="true to assume w_solute = w_solid at base of column, false to set w_solute=0.0 (ie zero solute velocity at great depth)"),
    )

end


function PB.register_methods!(rj::ReactionSedimentPhys)

    oceanfloor_vars = [  
        PB.VarDepColumn("oceanfloor_temp"=>"oceanfloor.temp",               "K",        "oceanfloor temperature"),
        PB.VarDepColumn("oceanfloor_sal"=>"oceanfloor.sal",                 "psu",      "oceanfloor salinity"),
        PB.VarDepColumn("oceanfloor_phi"=>"oceanfloor.phi",                 "",         "sediment surface porosity"),
        PB.VarDepColumn("oceanfloor_phimin"=>"(oceanfloor.phimin)",         "",         "sediment porosity at infinite depth"),
        PB.VarDepColumn("oceanfloor_w_accum"=>"(oceanfloor.w_accum)",       "m yr-1",   "sediment accumulation rate (+ve)"),
    ]

    phys_vars = [
        PB.VarProp("phi",                           "",         "cell mean porosity (volume fraction of solute phase)"),
        PB.VarProp("phi_solid",                     "",         "1.0 - cell mean porosity (volume fraction of solid phase)"),
        PB.VarProp("phi_upper",                     "",         "porosity (volume fraction of solute phase) at cell upper face"),
        PB.VarProp("phi_solid_upper",               "",         "1.0 - porosity (volume fraction of solid phase) at cell upper face"),
        PB.VarProp("phi_lower",                     "",         "porosity (volume fraction of solute phase) at cell lower face"),
        PB.VarProp("phi_solid_lower",               "",         "1.0 - porosity (volume fraction of solid phase) at cell lower face"),
        PB.VarProp("volume_solute",                 "m^3",      "solute volume of sediment cells"),
        PB.VarProp("volume_solid",                  "m^3",      "solid volume of sediment cells"),

        PB.VarProp("temp",                          "Kelvin",   "sediment temperature"),
        PB.VarProp("sal",                           "psu",      "sediment salinity"),
        PB.VarProp("Dfac",                          "",         "tortuoisity-dependent multiplier for solute diffusivity"),
    ]

    w_vars = [
        PB.VarTarget("volume_change_sms",           "m^3 yr-1", "change in volume due to production and consumption of solid phase components"),
        PB.VarProp("w_solid_upper",                 "m yr-1",   "solid phase advection velocity (downwards is -ve) across upper cell face"),
        PB.VarProp("w_solid_lower",                 "m yr-1",   "solid phase advection velocity (downwards is -ve) across lower cell face"),
        PB.VarProp("w_solute_upper",                "m yr-1",   "solute phase advection velocity (downwards is -ve) across upper cell face"),
        PB.VarProp("w_solute_lower",                "m yr-1",   "solute phase advection velocity (downwards is -ve) across lower cell face"),
    ]

    PB.add_method_do!(
        rj,
        do_sediment_phys,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(oceanfloor_vars),
            PB.VarList_namedtuple(phys_vars),
        ),
    )

    PB.add_method_do!(
        rj,
        do_sediment_w,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(oceanfloor_vars),
            PB.VarList_namedtuple(PB.VarDep.(phys_vars)),
            PB.VarList_namedtuple(w_vars),
        ),
    )

    # add to setup as well, so volumes are available for Variable initialization
    PB.add_method_setup!(
        rj,
        setup_sediment_phys,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(oceanfloor_vars),
            PB.VarList_namedtuple(phys_vars),
        ),
    )

    PB.add_method_initialize_zero_vars_default!(rj) # zero out volume_change_sms at start of each timestep

    return nothing
end


function setup_sediment_phys(
    m::PB.ReactionMethod,
    pars,
    vars,
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    attribute_name == :setup || return
    do_sediment_phys(m, pars, vars, cellrange, 0.0)
end

function do_sediment_phys(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, oceanfloor_vars, phys_vars),
    cellrange::PB.AbstractCellRange,
    deltat,
)

    for (icol, colindices) in cellrange.columns
        # porosity
        if pars.f_porosity[] == "Const"
            phys_vars.phi[colindices] .= oceanfloor_vars.oceanfloor_phi[icol]
            phys_vars.phi_upper[colindices] .= oceanfloor_vars.oceanfloor_phi[icol]
            phys_vars.phi_lower[colindices] .= oceanfloor_vars.oceanfloor_phi[icol]
        elseif pars.f_porosity[] == "ExpAtten"
            # exponential decline from surface to depth with lengthscale zpor
            phys_vars.phi[colindices] .= ((oceanfloor_vars.oceanfloor_phimin[icol] .+
                (oceanfloor_vars.oceanfloor_phi[icol] - oceanfloor_vars.oceanfloor_phimin[icol])
                .*exp.(@view(grid_vars.zmid[colindices])./pars.zpor[])))
            phys_vars.phi_upper[colindices] .= ((oceanfloor_vars.oceanfloor_phimin[icol] .+
                (oceanfloor_vars.oceanfloor_phi[icol] - oceanfloor_vars.oceanfloor_phimin[icol])
                .*exp.(@view(grid_vars.zupper[colindices])./pars.zpor[])))
            phys_vars.phi_lower[colindices] .= ((oceanfloor_vars.oceanfloor_phimin[icol] .+
                (oceanfloor_vars.oceanfloor_phi[icol] - oceanfloor_vars.oceanfloor_phimin[icol])
                .*exp.(@view(grid_vars.zlower[colindices])./pars.zpor[])))
        else
            error("unknown f_porosity='$(pars.f_porosity[])'")
        end

        for i in colindices
            phys_vars.temp[i]  =  oceanfloor_vars.oceanfloor_temp[icol]
            phys_vars.sal[i]   =  oceanfloor_vars.oceanfloor_sal[icol]

            phys_vars.phi_solid[i] = 1.0 - phys_vars.phi[i]
            phys_vars.phi_solid_upper[i] = 1.0 - phys_vars.phi_upper[i]
            phys_vars.phi_solid_lower[i] = 1.0 - phys_vars.phi_lower[i]
            phys_vars.volume_solute[i] = grid_vars.volume[i] * phys_vars.phi[i]
            phys_vars.volume_solid[i]  = grid_vars.volume[i] * phys_vars.phi_solid[i]

            # tortuoisity-dependent multiplier for diffusivity
            phys_vars.Dfac[i]  = 1.0 /(1.0 - log(phys_vars.phi[i]^2))  # Boundreau (1996) formulation

        end

    end

    return nothing
end


function do_sediment_w(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, oceanfloor_vars, phys_vars, w_vars),
    cellrange::PB.AbstractCellRange,
    deltat,
)

    for (icol, colindices) in cellrange.columns
        # indices for upper, lower cells 
        iu, il = first(colindices), last(colindices)

        # solid advection due to change in solid-phase volume
        last_w_solid_upper = zero(w_vars.w_solid_upper[iu])
        for i in colindices            
            w_vars.w_solid_upper[i] = last_w_solid_upper            
            # m yr-1 =   m3 yr-1          / m^3
            dw = w_vars.volume_change_sms[i]/(phys_vars.phi_solid_lower[i]*grid_vars.Abox[i])
            w_vars.w_solid_lower[i] = w_vars.w_solid_upper[i] - dw
            last_w_solid_upper = w_vars.w_solid_lower[i]
        end

        # add prescribed solid advection velocity, calculated from accumulation rate at sediment surface and assuming time-independent prescribed porosity
        if !isnothing(oceanfloor_vars.oceanfloor_w_accum)
            for i in colindices  
                w_vars.w_solid_upper[i] += (-oceanfloor_vars.oceanfloor_w_accum[icol] *
                                    (1.0 - phys_vars.phi[iu]) / (1.0 - phys_vars.phi_upper[i]))
                w_vars.w_solid_lower[i] += (-oceanfloor_vars.oceanfloor_w_accum[icol] *
                                    (1.0 - phys_vars.phi[iu]) / (1.0 - phys_vars.phi_lower[i]))
            end
        end

        # solute advection, assuming time-independent prescribed porosity
        if pars.w_solute[]
            # solute advection velocity, assuming solute is comoving with solid at base of column
            phi_floor = phys_vars.phi_lower[il]
            w_solute_floor = w_vars.w_solid_lower[il]                
            w_vars.w_solute_upper[colindices] .= (w_solute_floor * phi_floor) ./ @view(phys_vars.phi_upper[colindices])
            w_vars.w_solute_lower[colindices] .= (w_solute_floor * phi_floor) ./ @view(phys_vars.phi_lower[colindices])
        else
            # no solute advection
            w_vars.w_solute_upper[colindices] .= 0.0 # assume zero flux at great depth
            w_vars.w_solute_lower[colindices] .= 0.0 # assume zero flux at great depth
        end

    end

    return nothing
end



end # module
