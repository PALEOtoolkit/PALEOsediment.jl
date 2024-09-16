module SedimentPhys

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOsediment

import Printf

import Infiltrator # julia debugger

"""
    ReactionSedimentPhys

Sediment physical parameters

The sediment grid is defined by eg [`PALEOsediment.Sediment.SedimentGrid.ReactionSedimentGridn1D`](@ref).

## Physical environment

Accumulation rate, and porosity are set by `oceanfloor` Variables.


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


const oceanfloor_vars = [  
    PB.VarDepColumn("oceanfloor_temp"=>"oceanfloor.temp",               "K",        "oceanfloor temperature"),
    PB.VarDepColumn("oceanfloor_sal"=>"oceanfloor.sal",                 "psu",      "oceanfloor salinity"),
    PB.VarDepColumn("oceanfloor_phi"=>"oceanfloor.phi",                 "",         "sediment surface porosity"),
    PB.VarDepColumn("oceanfloor_phimin"=>"(oceanfloor.phimin)",         "",         "sediment porosity at infinite depth"),
    PB.VarDepColumn("oceanfloor_w_accum"=>"oceanfloor.w_accum",         "m yr-1",   "sediment accumulation rate (+ve)"),
]

const phys_vars = [
    PB.VarProp("phi",                           "",         "porosity (volume fraction of solute phase)"),
    PB.VarProp("phi_solid",                     "",         "1.0 - porosity (volume fraction of solid phase)"),
    PB.VarProp("volume_solute",                 "m^3",      "solute volume of sediment cells"),
    PB.VarProp("volume_solid",                  "m^3",      "solid volume of sediment cells"),

    PB.VarProp("temp",                          "Kelvin",   "sediment temperature"),
    PB.VarProp("sal",                           "psu",      "sediment salinity"),

    PB.VarProp("w_solid",                       "m yr-1",   "solid phase advection velocity (downwards is -ve)"),
    PB.VarProp("w_solute",                      "m yr-1",   "solute phase advection velocity (downwards is -ve)"),
    PB.VarProp("Dfac",                          "",         "tortuoisity-dependent multiplier for solute diffusivity"),
]

function PB.register_methods!(rj::ReactionSedimentPhys)

    PB.add_method_do!(
        rj,
        do_sediment_phys,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(oceanfloor_vars),
            PB.VarList_namedtuple(phys_vars),
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

"set porosity etc"
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
        elseif pars.f_porosity[] == "ExpAtten"
            # exponential decline from surface to depth with lengthscale zpor
            phys_vars.phi[colindices] .= ((oceanfloor_vars.oceanfloor_phimin[icol] .+
                (oceanfloor_vars.oceanfloor_phi[icol] - oceanfloor_vars.oceanfloor_phimin[icol])
                .*exp.(@view(grid_vars.zmid[colindices])./pars.zpor[])))
        else
            error("unknown f_porosity='$(pars.f_porosity[])'")
        end

        for i in colindices
            phys_vars.temp[i]  =  oceanfloor_vars.oceanfloor_temp[icol]
            phys_vars.sal[i]   =  oceanfloor_vars.oceanfloor_sal[icol]

            phys_vars.phi_solid[i] = 1.0 - phys_vars.phi[i]
            phys_vars.volume_solute[i] = grid_vars.volume[i] * phys_vars.phi[i]
            phys_vars.volume_solid[i]  = grid_vars.volume[i] * phys_vars.phi_solid[i]

            # solid advection velocity (solute advection calculated below)
            phys_vars.w_solid[i] = (-oceanfloor_vars.oceanfloor_w_accum[icol] *
                                (phys_vars.phi_solid[first(colindices)]) / (phys_vars.phi_solid[i]))
            
            ### NEW advection equation
           
            # uf (m y-1) displacement velocity (distance by which the cells are displaced per unit time)
            # of a cell at index i is equal to the added net specific mass production of all microbial species of the 
            # biofilm matrix between the substratum and this location

            # uf[i] = sum of R_sms (m3 yr-1) for index i 0-i divided by cell area A (m2)
            # where R_sms is PB.VarTarget("volume_change_sms",  "m^3 yr-1", "change in volume due to production and consumption of solid phase components")
            
            # w_solid adds sediment accumulation w_accum (m y-1) of solid phases to the solid phase biological accumulation uf,i
            # phys_vars.w_solid[i] = (-oceanfloor_vars.oceanfloor_w_accum[icol] *
            #                   (phys_vars.phi_solid[first(colindices)]) / (phys_vars.phi_solid[i]))
            #                   + uf[i]

            # tortuoisity-dependent multiplier for diffusivity
            phys_vars.Dfac[i]  = 1.0 /(1.0 - log(phys_vars.phi[i]^2))  # Boundreau (1996) formulation

        end

        if pars.w_solute[]
            # solute advection velocity, assuming solute is comoving with solid at base of column
            phi_floor = phys_vars.phi[last(colindices)]
            w_solute_floor = phys_vars.w_solid[last(colindices)]                
            phys_vars.w_solute[colindices] .= (w_solute_floor * phi_floor) ./ @view(phys_vars.phi[colindices])
        else
            # no solute advection
            phys_vars.w_solute[colindices] .= 0.0 # assume zero flux at great depth
        end

    end

    return nothing
end



end # module
