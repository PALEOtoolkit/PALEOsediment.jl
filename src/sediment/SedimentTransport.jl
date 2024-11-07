module SedimentTransport

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOsediment
import PALEOaqchem

import LinearAlgebra
import SparseArrays
import Printf

import Infiltrator # julia debugger

#########################################################################################################
# Sediment solid phase transport
#########################################################################################################

"""
    ReactionSedimentTransportSolid

Sediment solid-phase transport for n x 1D sediment columns.

The sediment grid is defined by eg [`PALEOsediment.Sediment.SedimentGrid.ReactionSedimentGridn1D`](@ref).

Physical parameters (porosity, advection velocity) are defined by eg [`PALEOsediment.Sediment.SedimentPhys.ReactionSedimentPhys](@ref).

## Physical environment

Bioturbation rates should be supplied on the sediment grid by eg [`PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates`](@ref).

## Transport components and species

Each component `<totalname>` to be transported should be defined by a source-minus-sink flux `<totalname>_sms`,
and one or more concentration Variables with names of form `<speciesname>_conc`

Solid concentration Variables are identified by attributes `vphase == VP_Solid` and `advect == true`,
and are transported by bioturbation and advection.

Transport fluxes are then accumulated into vectors `<totalname>_sms`:
- If the `totalnames` attribute is not set  on a variable `<speciesname>_conc`, then `totalname` is assumed to be the 
  same as `speciesname` and transport fluxes are accumulated into `<speciesname>_sms`.
- If the `totalnames` attribute is set to a Vector of strings `["mult_1*totalname_1", ... "mult_n*totalname_n"]`, 
  then this is used to define the appropriate vector of `[<totalname_n>_sms]` fluxes (allowing multiple species 
  concentrations with different transport properties or phases for each `<totalname>_sms`, and each species to contribute 
  to multiple `<totalname_n>_sms` with a stoichiometry multiplier `mult_n`).

## Boundary conditions

Input particulate fluxes should be added by the `fluxOceanfloor` flux coupler to the surface sediment cells
by linking to Variables in the `sediment.oceanfloor` subdomain, ie to `sediment.oceanfloor.<totalname>_sms`.

Burial fluxes are the base of the sediment columns are defined in the `Domain` `fluxOceanBurial`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedimentTransportSolid{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        # PB.ParBool("chargebalance", false,
        #    description="TODO true to calculate Coulombic coupling and migration to maintain charge balance"),

    )


    # solid_charges::Vector{Float64} = Float64[] # TODO
    solid_totals_multipliers::Vector{Vector{Float64}} = Vector{Float64}[]
end


function PB.register_methods!(rj::ReactionSedimentTransportSolid)
    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionSedimentTransportSolid)

    phys_vars = [
        PB.VarDep("phi_solid",                     "",         "cell mean 1.0 - porosity (volume fraction of solid phase)"),
        PB.VarDep("phi_solid_upper",               "",         "1.0 - porosity (volume fraction of solid phase) on cell upper face"),
        PB.VarDep("phi_solid_lower",               "",         "1.0 - porosity (volume fraction of solid phase) on cell lower face"),
        PB.VarDep("volume_solid",                  "m^3",      "solid volume of sediment cells"),
        PB.VarDep("w_solid_upper",                 "m yr-1",   "solid phase advection velocity (downwards is -ve) at cell upper face"),
        PB.VarDep("w_solid_lower",                 "m yr-1",   "solid phase advection velocity (downwards is -ve) at cell lower face"),
    ]


    bio_vars = [
        PB.VarDep("diff_bioturb",                 "m^2 yr-1", "bioturbation effective diffusivity"),
    ]

    (
        vars_solid_conc,
        vvars_solid_sms,
        vvars_solid_fluxOceanBurial,
        totals_mult,
    ) = find_solid_vars(rj.domain)
    empty!(rj.solid_totals_multipliers)
    append!(rj.solid_totals_multipliers, totals_mult)

    PB.add_method_do!(
        rj,
        do_sediment_transport_solid,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(phys_vars),
            PB.VarList_namedtuple(bio_vars),
            PB.VarList_tuple(vars_solid_conc),
            PB.VarList_ttuple(vvars_solid_sms),
            PB.VarList_ttuple(vvars_solid_fluxOceanBurial),
        ),
    )

    # output transport conc and flux variable list to log
    PB.add_method_setup!(
        rj,
        setup_sediment_transport_solid,
        (
            PB.VarList_vector(vars_solid_conc),
            PB.VarList_vvector(vvars_solid_sms),
        )
    )

    return nothing
end


function find_solid_vars(domain::PB.Domain)
    # Find all solid phase concentration variables with attribute :vphase == VP_Solid, and :advect == true,
    # and name of form `<rootname>_conc`.
    # Define `totalname` from :totalnames attribute if present, otherwise set `totalname = rootname`.
    # Then provide links to `<totalname>_sms`, `fluxOceanBurial.flux_<totalname>`.

    filter_conc(v) = (
        (PB.get_attribute(v, :vphase, PB.VP_Undefined) == PB.VP_Solid)
        && PB.get_attribute(v, :advect, false)
    )
    conc_domvars = PB.get_variables(domain, filter_conc)

    vars_conc, vvars_sms, vvars_fluxOceanBurial, totals_mult = [], [], [], Vector{Float64}[]
    for v in conc_domvars
        v.name[end-4:end] == "_conc" ||
            error("find_solid_vars: Variable $(PB.fullname(v)) has :advect attribute == true but is not named _conc")
        rootname = v.name[1:end-5]
       
        # per-species concentration
        # default :field_data to link to any ScalarData, IsotopeData etc
        push!(vars_conc,
            PB.VarDep(    rootname*"_conc", "mol m-3", "")
        )

        # per-totalvar _sms and flux output
        totalnamesmults = PB.get_attribute(v, :totalnames, missing)
        totalnamesmults = (ismissing(totalnamesmults) || isempty(totalnamesmults)) ? [rootname] : totalnamesmults
      
        # TODO look up and use existing variables, if present
        st_sms, st_fluxoceanburial, st_mult = [], [], Float64[] 
        for tnm in totalnamesmults
            tmult, tname = PALEOaqchem.parse_number_name(tnm)
            push!(st_sms, PB.VarContrib(tname*"_sms", "mol yr-1", ""))
            push!(st_fluxoceanburial,
                PB.VarContribColumn("fluxOceanBurial_flux_"*tname=>"fluxOceanBurial.flux_"*tname, "mol yr-1", "")
            )
            push!(st_mult, tmult)
        end
        push!(vvars_sms, st_sms)
        push!(vvars_fluxOceanBurial, st_fluxoceanburial)
        push!(totals_mult, st_mult)
    end

    return (vars_conc, vvars_sms, vvars_fluxOceanBurial, totals_mult)
end

function setup_sediment_transport_solid(
    m::PB.ReactionMethod,
    pars,
    (
        vars_solid_conc,
        vars_solid_sms,
    ),
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    attribute_name == :setup || return
    rj = m.reaction

    io = IOBuffer()
    println(io, "setup_sediment_transport_solid $(PB.fullname(rj)) ReactionSedimentTransportSolid")
    println(io)
    # Vector of VariableReactions corresponding to the (vars_solid_conc,) data arrays
    (rvars_solid_conc, rvars_solid_sms) = PB.get_variables_tuple(m; flatten=false)

    Printf.@printf(io, "    %30s%70s\n", "Solid Variable", "transport flux")
    for (rv_sc, rvs_sms, totals_mult) in PB.IteratorUtils.zipstrict(rvars_solid_conc, rvars_solid_sms, rj.solid_totals_multipliers)
        total_names_mult_str = string(["$m*$(PB.fullname(v.linkvar))" for (m, v) in PB.IteratorUtils.zipstrict(totals_mult, rvs_sms)]) 
        Printf.@printf(io, "    %30s%70s\n", PB.fullname(rv_sc.linkvar), total_names_mult_str)
    end

    @info String(take!(io))

    return nothing
end


function do_sediment_transport_solid(
    m::PB.ReactionMethod,
    (
        grid_vars,
        phys_vars,
        bio_vars,
        vars_solid_conc,
        vvars_solid_sms,
        vvars_solid_fluxOceanBurial,
    ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    function _calc_solid_species_transport(
        solid_conc,
        solid_totals_multipliers,
        solid_smss, # Tuple of arrays
        solid_fluxOceanBurials, # Tuple of arrays
        (
            grid_vars,
            phys_vars,
            bio_vars,
            icol, colindices,
        )
    )        
        calc_downwards_advect_column(
            grid_vars.Abox, phys_vars.phi_solid_upper, phys_vars.w_solid_upper, phys_vars.phi_solid_lower, phys_vars.w_solid_lower,
            solid_conc, solid_totals_multipliers, solid_smss,
            nothing, nothing, solid_fluxOceanBurials,
            icol, colindices
        )

        calc_diffuse_column(
            grid_vars.Abox, grid_vars.zupper, grid_vars.zmid, phys_vars.phi_solid,
            bio_vars.diff_bioturb, nothing, solid_conc, solid_totals_multipliers, solid_smss,
            colindices
        )

        return nothing
    end

    for (icol, colindices) in cellrange.columns
        PB.IteratorUtils.foreach_longtuple_p(
            _calc_solid_species_transport,
            vars_solid_conc,
            rj.solid_totals_multipliers,
            vvars_solid_sms,
            vvars_solid_fluxOceanBurial,
            (
                grid_vars,
                phys_vars,
                bio_vars,
                icol, colindices,
            ),
        )
    end

    return nothing
end

###############################################################################
# Sediment solute transport
###############################################################################

"""
    SedimentTransportSolute

Sediment solute transport for n x 1D sediment columns.

The sediment grid is defined by eg [`PALEOsediment.Sediment.SedimentGrid.ReactionSedimentGridn1D`](@ref).

Physical parameters (porosity, advection velocity) are defined by eg [`PALEOsediment.Sediment.SedimentPhys.ReactionSedimentPhys](@ref).

## Physical environment

Bioturbation and bioirrigation rates should be supplied on the sediment grid by eg [`PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates`](@ref).

## Transport components and species

Each component `<totalname>` to be transported should be defined by a source-minus-sink flux `<totalname>_sms`,
and one or more concentration Variables with names of form `<speciesname>_conc`

Solute concentration Variables are identified by attributes `vphase == VP_Solute` and `advect == true`,
and are transported by diffusion, bioturbation and bioirrigation, and advection.
Species-specific solute diffusivities are calculated either based on the attribute `diffusivity_speciesname` of the
`<speciesname>_conc` solute Variables, which should be one of the  names available from
`PALEOaqchem.MolecularDiffusion.create_solute_diffusivity_func`, or (if `diffusivity_speciesname` is not set ie an empty string)
by the attribute `diffusivity` (units cm^2 s-1).

Transport fluxes are then accumulated into vectors `<totalname>_sms`:
- If the `totalnames` attribute is not set  on a variable `<speciesname>_conc`, then `totalname` is assumed to be the 
  same as `speciesname` and transport fluxes are accumulated into `<speciesname>_sms`.
- If the `totalnames` attribute is set to a Vector of strings `["mult_1*totalname_1", ... "mult_n*totalname_n"]`, 
  then this is used to define the appropriate vector of `[<totalname_n>_sms]` fluxes (allowing multiple species 
  concentrations with different transport properties or phases for each `<totalname>_sms`, and each species to contribute 
  to multiple `<totalname_n>_sms` with a stoichiometry multiplier `mult_n`).

## Boundary conditions

Oceanfloor solute fluxes should be defined in the `Domain` `fluxOceanfloor`, with names `fluxOceanfloor.soluteflux_<totalname>`.

Burial fluxes are the base of the sediment columns are defined in the `Domain` `fluxOceanBurial`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedimentTransportSolute{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParDouble("zdbl", 0.04e-2, units="m",
            description="diffusive boundary layer thickness at sediment - water column interface"),

        PB.ParBool("solute_burial_flux", false,
            description="true to include solute burial flux terms"),

        # PB.ParBool("chargebalance", false,
        #    description="TODO true to calculate Coulombic coupling and migration to maintain charge balance"),

    )


    solute_gammas::Vector{Float64} = Float64[]
    # solute_charges::Vector{Float64} = Float64[] # TODO
    solute_totals_multipliers::Vector{Vector{Float64}} = Vector{Float64}[]
end



function PB.register_methods!(rj::ReactionSedimentTransportSolute)
    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionSedimentTransportSolute)

    phys_vars = [
        PB.VarDep("phi",                           "",         "cell mean porosity (volume fraction of solute phase)"),
        PB.VarDep("volume_solute",                 "m^3",      "solute volume of sediment cells"),
        PB.VarDep("phi_upper",                     "",         "cell mean porosity (volume fraction of solute phase) at cell upper face"),
        PB.VarDep("phi_lower",                     "",         "cell mean porosity (volume fraction of solute phase) at cell lower face"),

        PB.VarDep("temp",                          "Kelvin",   "sediment temperature"),
        PB.VarDep("sal",                           "psu",      "sediment salinity"),

        PB.VarDep("w_solute_upper",                "m yr-1",   "solute phase advection velocity (downwards is -ve) at cell upper face"),
        PB.VarDep("w_solute_lower",                "m yr-1",   "solute phase advection velocity (downwards is -ve) at cell lower face"),
        PB.VarDep("Dfac",                          "",         "tortuoisity-dependent multiplier for solute diffusivity"),
    ]


    bio_vars = [
        PB.VarDep("diff_bioturb",                 "m^2 yr-1", "bioturbation effective diffusivity"),
        PB.VarDep("alpha_bioirrig",               "yr-1",     "bioirrigation exchange rate"),
    ]


    (
        vars_solute_conc,
        vvars_solute_sms,
        vars_solute_oceanfloor_conc,
        vvars_solute_fluxOceanfloor,
        vvars_solute_fluxOceanBurial,  # vector of empty vectors if rj.pars_w_solute[] == false
        vars_solute_diff,
        totals_mult,
    ) = find_solute_vars(rj.domain; include_burial_flux=rj.pars.solute_burial_flux[])
    empty!(rj.solute_totals_multipliers)
    append!(rj.solute_totals_multipliers, totals_mult)

    PB.add_method_do!(
        rj,
        do_sediment_transport_solute,
        (
            PB.VarList_namedtuple(PB.VarDep.(PALEOsediment.Sediment.SedimentGrid.grid_vars)),
            PB.VarList_namedtuple(phys_vars),
            PB.VarList_namedtuple(bio_vars),
            PB.VarList_tuple(vars_solute_conc),
            PB.VarList_ttuple(vvars_solute_sms),
            PB.VarList_tuple(vars_solute_oceanfloor_conc),
            PB.VarList_ttuple(vvars_solute_fluxOceanfloor),
            PB.VarList_ttuple(vvars_solute_fluxOceanBurial),
            PB.VarList_tuple(vars_solute_diff),
            # fns_solutediff added in preparefn
        ),
        preparefn=prepare_do_sediment_transport_solute,
    )


    # read solute variable attributes, output transport conc and flux variable list to log
    PB.add_method_setup!(
        rj,
        setup_sediment_transport_solute,
        (
            PB.VarList_vector(vars_solute_conc),
            PB.VarList_vvector(vvars_solute_sms),
        )
    )

    return nothing
end

"""
    find_solute_vars(domain::Domain; include_burial_flux) -> (vars_conc, vvars_sms, vars_oceanfloor_conc, vvars_fluxOceanfloor, vars_diff, totals_mult)

Find all solute concentration variables with attribute :vphase == VP_Solute, and :advect == true,
and name of form `<rootname>_conc`.
Define `totalnames` from :totalnames attribute if present, otherwise set `totalnames = [rootname]`.
Then provide links to `[<totalnames>_sms]`, `oceanfloor.<rootname>_conc`, `[fluxOceanfloor.soluteflux_<totalnames>]`
"""
function find_solute_vars(domain::PB.Domain; include_burial_flux)

    filter_conc(v) = (
        (PB.get_attribute(v, :vphase, PB.VP_Undefined) == PB.VP_Solute)
        && PB.get_attribute(v, :advect, false)
    )
    conc_domvars = PB.get_variables(domain, filter_conc)

    vars_conc, vvars_sms, vars_oceanfloor_conc, vvars_fluxOceanfloor, vvars_fluxOceanBurial, vars_diff, totals_mult = [], [], [], [], [], [], Vector{Float64}[]
    for v in conc_domvars
        v.name[end-4:end] == "_conc" ||
            error("find_solute_vars: Variable $(PB.fullname(v)) has :advect attribute == true but is not named _conc")
        rootname = v.name[1:end-5]
        
        # per-species sediment concentration and upper boundary condition

        # default :field_data to link to any ScalarData, IsotopeData etc
        push!(vars_conc,
            PB.VarDep(    rootname*"_conc", "mol m-3", "")
        )
        
        push!(vars_oceanfloor_conc,
            PB.VarDepColumn(    "oceanfloor_"*rootname*"_conc"=>"oceanfloor."*rootname*"_conc", "mol m-3", "")
        )

        push!(vars_diff,
            PB.VarProp(   rootname*"_diff", "m^2/yr", "solute diffusivity")
        )

        # per-totalvar _sms and flux output
        totalnamesmults = PB.get_attribute(v, :totalnames, missing)
        totalnamesmults = (ismissing(totalnamesmults) || isempty(totalnamesmults)) ? [rootname] : totalnamesmults
      
        # TODO look up and use existing variable, if present
        st_sms, st_fluxoceanfloorsolute, st_fluxoceanburial, st_mult = [], [], [], Float64[] 
        for tnm in totalnamesmults
            tmult, tname = PALEOaqchem.parse_number_name(tnm)
            push!(st_sms, PB.VarContrib(tname*"_sms", "mol yr-1", ""))
            push!(st_fluxoceanfloorsolute, PB.VarContribColumn("fluxOceanfloor_soluteflux_"*tname=>"fluxOceanfloor.soluteflux_"*tname, "mol yr-1", ""))
            if include_burial_flux
                push!(st_fluxoceanburial, PB.VarContribColumn("fluxOceanBurial_flux_"*tname=>"fluxOceanBurial.flux_"*tname, "mol yr-1", ""))
            end
            push!(st_mult, tmult)
        end
        push!(vvars_sms, st_sms)
        push!(vvars_fluxOceanfloor, st_fluxoceanfloorsolute)
        push!(vvars_fluxOceanBurial, st_fluxoceanburial)
        push!(totals_mult, st_mult)

    end

    return (vars_conc, vvars_sms, vars_oceanfloor_conc, vvars_fluxOceanfloor, vvars_fluxOceanBurial, vars_diff, totals_mult)
end


"read :diffusivity_speciesname or :diffusivity attribute from _conc Variable and create solute diffusivity functions"
function prepare_do_sediment_transport_solute(m::PB.ReactionMethod, vardata)

    (_, _, _, vars_solute_conc, _, _, _, _, _,) = PB.get_variables_tuple(m)

    diff_fns = []
    for v_conc in vars_solute_conc
        # solute diffusivity accessors and calculation method
        diffusivity_speciesname = PB.get_domvar_attribute(v_conc, :diffusivity_speciesname, missing)
        if ismissing(diffusivity_speciesname) || isempty(diffusivity_speciesname)
            diffusivity = PB.get_domvar_attribute(v_conc, :diffusivity, missing)
            if diffusivity isa Float64 
                diffusivity_speciesname = string(diffusivity)
            end
        end
        if ismissing(diffusivity_speciesname) || isempty(diffusivity_speciesname)
            @error("prepare_do_sediment_transport_solute: no :diffusivity_speciesname or :diffusivity attribute found for Variable $(PB.fullname(v_conc.linkvar))")
        end

        push!(diff_fns, PALEOaqchem.MolecularDiffusion.create_solute_diffusivity_func(diffusivity_speciesname))
    end

    return (vardata..., Tuple(diff_fns))
end

function setup_sediment_transport_solute(
    m::PB.ReactionMethod,
    pars,
    (
        vars_solute_conc,
        vars_solute_sms,
    ),
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    attribute_name == :setup || return
    rj = m.reaction

    io = IOBuffer()
    println(io, "setup_sediment_transport $(PB.fullname(rj)) ReactionSedimentTransportSolute")
    println(io)
    # Vector of VariableReactions corresponding to the (vars_solute_conc,) data arrays
    (rvars_solute_conc, rvars_solute_sms) = PB.get_variables_tuple(m; flatten=false)

    Printf.@printf(io, "    %30s%20s%50s\n", "Solute Variable", "gamma", "transport flux")
    empty!(rj.solute_gammas)
    for (rv_sc, rvs_sms, totals_mult) in PB.IteratorUtils.zipstrict(rvars_solute_conc, rvars_solute_sms, rj.solute_totals_multipliers)
        conc_domvar = rv_sc.linkvar # read attribute from linked Domain variable 
        gamma = coalesce(PB.get_attribute(conc_domvar, :gamma, missing), 1.0) # either attribute not present, or present with value missing -> 1.0        
        total_names_mult_str = string(["$m*$(PB.fullname(v.linkvar))" for (m, v) in PB.IteratorUtils.zipstrict(totals_mult, rvs_sms)])
        Printf.@printf(io, "    %30s%20g%50s\n", PB.fullname(conc_domvar), gamma, total_names_mult_str)
        if !isa(gamma, Float64)
            @info String(take!(io))
            error("read invalid value $gamma ($(typeof(gamma)) from $(PB.fullname(conc_domvar))")
        end
        push!(rj.solute_gammas, gamma)
    end

    @info String(take!(io))

    return nothing
end

function do_sediment_transport_solute(
    m::PB.ReactionMethod,
    pars,
    (
        grid_vars,
        phys_vars,
        bio_vars,
        vars_solute_conc,
        vvars_solute_sms,
        vars_solute_oceanfloor_conc,
        vvars_solute_fluxOceanfloor,
        vvars_solute_fluxOceanBurial,
        vars_solute_diff,
        fns_solutediff, # added by prepare_
    ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    # @Infiltrator.infiltrate

    function _calc_solute_species_transport(
        solute_conc,
        solute_totals_multipliers, 
        solute_smss,  # Tuple of arrays
        solute_oceanfloor_conc,
        solute_fluxOceanfloors, # Tuple of arrays
        solute_fluxOceanBurial, # Tuple of arrays
        solute_diff,
        solute_gamma,
        (
            zdbl,
            grid_vars,
            phys_vars,
            bio_vars,
            icol, colindices,
        ),
    )    
        calc_downwards_advect_column(
            grid_vars.Abox, phys_vars.phi_upper, phys_vars.w_solute_upper, phys_vars.phi_lower, phys_vars.w_solute_lower,
            solute_conc, solute_totals_multipliers, solute_smss,
            solute_oceanfloor_conc,  solute_fluxOceanfloors, solute_fluxOceanBurial,
            icol, colindices
        )
    
        calc_diffuse_column(
            grid_vars.Abox, grid_vars.zupper, grid_vars.zmid, phys_vars.phi,
            solute_diff, bio_vars.diff_bioturb, solute_conc, solute_totals_multipliers, solute_smss,
            colindices
        )
    
        calc_surface_dbl_column(
            grid_vars.Abox, zdbl, phys_vars.phi_upper,
            solute_diff, phys_vars.Dfac, solute_conc, solute_totals_multipliers, solute_smss,
            solute_oceanfloor_conc, solute_fluxOceanfloors,
            icol, colindices
        )
    
        calc_irrigate_column(
            phys_vars.volume_solute,
            bio_vars.alpha_bioirrig, solute_gamma, solute_conc, solute_totals_multipliers, solute_smss,
            solute_oceanfloor_conc, solute_fluxOceanfloors,
            icol, colindices
        )
    
    
        return nothing
    end

    for (icol, colindices) in cellrange.columns
        PB.IteratorUtils.foreach_longtuple_p(
            calc_solute_species_diffusivity,
            vars_solute_diff,
            fns_solutediff,
            (
                phys_vars.temp,
                grid_vars.pressure,
                phys_vars.sal,
                phys_vars.Dfac,
                colindices,
            ),
        )

        PB.IteratorUtils.foreach_longtuple_p(
            _calc_solute_species_transport,
            vars_solute_conc,
            rj.solute_totals_multipliers,
            vvars_solute_sms,
            vars_solute_oceanfloor_conc,
            vvars_solute_fluxOceanfloor,
            vvars_solute_fluxOceanBurial,
            vars_solute_diff,
            rj.solute_gammas,
            (
                pars.zdbl[],
                grid_vars,
                phys_vars,
                bio_vars,
                icol, colindices,
            ),
        )

    end

    return nothing
end


"calculate tortuoisity-corrected solute diffusivity"
function calc_solute_species_diffusivity(
    solute_diff,
    f_solute_diff,
    (
        temp,
        pressure,
        sal,
        Dfac,
        indices
    )
)

    @inbounds for i in indices
        #   m^2 yr-1             cm^2 s-1      bar dbar-1 dbar                      s yr-1  (cm m-1)^2
        solute_diff[i] = (
            f_solute_diff(temp[i], 0.1*pressure[i]+1.0, sal[i])*PB.Constants.k_secpyr*1e-4*Dfac[i]
        )
    end
    return nothing
end

"calculate sediment surface flux across diffusive boundary layer thicknes zdbl"
function calc_surface_dbl_column(
    area, zdbl, volfrac,
    diff_tortuoisity, Dfac, conc, mults, smss,
    ub_conc, ub_fluxouts,
    icol, colindices
)

    i = first(colindices)
    # exchange at sediment-water interface, with diffusive boundary layer thickness zdbl
    dz = zdbl           # assume diffusion across zdbl to a sediment top cell with uniform concentration
    dc = ub_conc[icol] - conc[i]
    sigma = diff_tortuoisity[i]/Dfac[i] # undo tortuosity correction to solute diffusivity
    a = area[i]*volfrac[i]  # area * porosity
    flux = -a*sigma*dc/dz
    # ub_fluxout[icol] += mult*flux
    _add_fluxes(ub_fluxouts, mults, (icol, flux))
    # sms[i] -= mult*flux
    _add_fluxes(smss, mults, (i, -flux))

    return nothing
end

"calculate bioirrigation fluxes for a single column"
function calc_irrigate_column(
    vol,
    alpha, solute_gamma, conc, mults, smss,
    ub_conc, ub_fluxouts,
    icol, colindices
)

    for i in colindices
    # mol yr-1   yr-1   mol m-3                 m^3
        flux = solute_gamma*alpha[i]*(ub_conc[icol] - conc[i])*vol[i]
       
        # sms[i] += mult*flux
        _add_fluxes(smss, mults, (i, flux))
        # ub_fluxout[icol] -= mult*flux
        _add_fluxes(ub_fluxouts, mults, (icol, -flux))
    end

    return nothing
end



#############################################################################################################
# Sediment column transport helper functions
#############################################################################################################

"downwards advection in a column at velocity  w (-ve)

 icol is column index for ub_conc, lb_fluxouts. colindices is column indices (top -> bottom)
for area, volfrac, w, conc, smss.

If `ub_conc = nothing`, no flux into upper boundary.

If `lb_fluxouts = nothing`, no flux out of lower boundary, flux accumulates in lowest cell"
function calc_downwards_advect_column(
    area, volfrac_upper, w_upper, volfrac_lower, w_lower,
    conc, mults, smss,
    ub_conc, ub_fluxins, lb_fluxouts,
    icol, colindices
)
    iu = first(colindices)
    if isnothing(ub_conc)
        flux = zero(conc[iu])
    else
        # flux (+ve) downwards into cell at top of column
        # NB: add zero(conc[iu]) to maintain type stability 
        # mol yr-1             mol m-3                 m yr-1   m^2
        flux = zero(conc[iu]) + ub_conc[icol]*volfrac_upper[iu]*(-w_upper[iu])*area[iu]       
        _add_fluxes(ub_fluxins, mults, (icol, -flux)) # +ve is out of sediment
    end
    @inbounds for i in colindices
        # add flux from above to this cell        
        _add_fluxes(smss, mults, (i, flux))

        # flux (+ve) leaving lower boundary of cell
        w_lower[i] <= 0.0 || error("upwards advection not supported")
        # mol yr-1             mol m-3                 m yr-1   m^2
        flux                = conc[i] * volfrac_lower[i] * (-w_lower[i]) * area[i]
        _add_fluxes(smss, mults, (i, -flux))
    end

    if isnothing(lb_fluxouts) || isempty(lb_fluxouts)
        # no flux out of lower boundary, add back to lowest cell
        _add_fluxes(smss, mults, (last(colindices), flux))
    else
        _add_fluxes(lb_fluxouts, mults, (icol, flux))
    end

    return nothing
end

"diffusive transport in a column interior with no-flux boundary conditions

colindices are column indices (top -> bottom)
for area, zupper, zmid, volfrac, two contributions to diff, conc, smss.

TODO - use phi at cell face instead of interpolating from mid-points ?"
function calc_diffuse_column(
    area, zupper, zmid, volfrac,
    diff1, diff2, conc, mults, smss,
    colindices,
)

    lasti = zero(first(colindices)) # set type of lasti
    @inbounds for i in colindices # count down in column from top
        if i != first(colindices)
            # diffusion in interior, with no-flux boundary at base of column
            dz = zmid[lasti] - zmid[i]
            dc = conc[lasti] - conc[i]
            # linearly interpolate diff, area
            wtu = (zmid[lasti] - zupper[i])/dz
            wtl = (zupper[i]-zmid[i])/dz
            sigma = wtu*diff1[lasti] + wtl*diff1[i]
            if !isnothing(diff2)
                sigma += wtu*diff2[lasti] + wtl*diff2[i]
            end
            a = wtu*area[lasti]*volfrac[lasti] + wtl*area[i]*volfrac[i]
            flux = -a*sigma*dc/dz # +ve is upwards in column
            _add_fluxes(smss, mults, (lasti, flux))
            _add_fluxes(smss, mults, (i, -flux))
        end

        lasti = i
    end

    return nothing
end



"""
    _add_fluxes(vecs, mults, (idx, flux))

Add flux * mult to element idx of each of a Tuple of vecs, equivalent to:    

    for k in 1:length(vecs)
        vecs[k][idx] += mults[k]*flux
    end
"""
@inline _add_fluxes(vecs, mults, (idx, flux)) = PB.IteratorUtils.foreach_tuple_unchecked_p(_add_flux, vecs, mults, (idx, flux))
@inline _add_flux(vec, m, (idx, flux))  = vec[idx] += m*flux


end # module
