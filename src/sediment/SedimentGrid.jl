module SedimentGrid

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOsediment

import Printf

import Infiltrator # julia debugger

const grid_vars = [
    PB.VarPropStateIndep("volume",              "m^3",      "volume of sediment cells"),
    PB.VarPropScalarStateIndep("volume_total",  "m^3",      "total volume of sediment cells"),
    PB.VarPropStateIndep("Abox",                "m^2",      "horizontal area of box"),
    PB.VarPropStateIndep("zupper",              "m",        "depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m"),
    PB.VarPropStateIndep("zlower",              "m",        "depth of lower surface of box (m)"),
    PB.VarPropStateIndep("zmid",                "m",        "mean depth of box"),
    PB.VarPropStateIndep("pressure",            "dbar",     "sediment pressure"),
    PB.VarPropStateIndep("rho_ref",             "kg m^-3",  "density conversion factor"),
]

const oceanfloor_grid_vars = [
    PB.VarDepColumnStateIndep("oceanfloor_Afloor"=>"oceanfloor.Afloor", "m^2",      "horizontal area of seafloor at sediment surface"),
    PB.VarDepColumnStateIndep("oceanfloor_zfloor"=>"oceanfloor.zfloor", "m",        "depth of ocean floor (m, -ve)"),
]


"""
    ReactionSedimentGridn1D

Sediment grid for n x 1D sediment columns.


A grid with n columns is created in a `Domain` `sediment`, bounded at the top by `Domain` `oceanfloor` and
at the base by `Domain` `sedimentfloor`.  Boundary cells at the sediment surface are therefore
in subdomain `sediment.oceanfloor`.

The number of columns, and column area are set by `oceanfloor` Variables.


# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedimentGridn1D{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("L",   0.15,   units="m",
            description="depth of sediment column"),
        PB.ParInt("ncellspercol",  60,
            description="number of cells per column"),
        PB.ParString("f_grid", "linear", allowed_values=["linear", "quadratic"],
            description="vertical grid transformation"),
        PB.ParDouble("grid_eta", NaN, units="m",
            description="length scale for vertical grid transformation"),

        PB.ParDouble("rho_ref", 1027.0, units="kg m-3",
            description = "assumed constant sw density conversion factor"),
    )


    "boundary and fluxTransfer Domains"
    domain_oceanfloor::Union{PB.Domain, Nothing}       = nothing
    domain_fluxOceanfloor::Union{PB.Domain, Nothing}   = nothing
    domain_sedimentfloor::Union{PB.Domain, Nothing}    = nothing
    domain_fluxOceanBurial::Union{PB.Domain, Nothing}  = nothing

end

# Define sediment grid for n x 1D columns, where n is derived from Oceanfloor Domain
function PB.set_model_geometry(rj::ReactionSedimentGridn1D, model::PB.Model)

    rj.domain_oceanfloor = PB.get_domain(model, "oceanfloor")
    !isnothing(rj.domain_oceanfloor) ||
        error("SedimentTransport.set_model_geometry no oceanfloor Domain")

    # size may not have been set yet
    !isnothing(rj.domain_oceanfloor.grid) ||
        error("Sediment.transport.set_model_geometry oceanfloor no grid set")

    ncolumns = rj.domain_oceanfloor.grid.ncells
    ncellspercol = rj.pars.ncellspercol[]
    @info "ReactionSedimentGridn1D set_model_geometry: $(PB.fullname(rj))"
    @info "    Defining $ncolumns sediment columns each $ncellspercol cells"

    # create sediment grid
    rj.domain.grid = PB.Grids.UnstructuredColumnGrid(
        ncells=ncolumns*ncellspercol,
        Icolumns=[collect((ic-1)*ncellspercol+1:ic*ncellspercol) for ic in 1:ncolumns]
    )

    # set subdomain mappings
    isedsurf = [col[1] for col in rj.domain.grid.Icolumns]  # first index in each column is surface box
    isedfloor = [col[end] for col in rj.domain.grid.Icolumns]  # last index in each column is surface box

    PB.Grids.set_subdomain!(rj.domain.grid, "oceanfloor", PB.Grids.BoundarySubdomain(isedsurf), true)
    @info "  set sediment.sedimentsurface Subdomain size=$(length(isedsurf))"
    PB.Grids.set_subdomain!(rj.domain.grid, "sedimentfloor", PB.Grids.BoundarySubdomain(isedfloor), true)
    @info "  set sediment.sedimentfloor Subdomain size=$(length(isedfloor))"

    # set sediment.grid
    domain_sediment = PB.get_domain(model, "sediment")
    !isnothing(domain_sediment) ||
        error("SedimentTransport.set_model_geometry no sediment Domain")
    domain_sediment.grid = rj.domain.grid

    # if not already set by Ocean component, set oceanfloor.grid, fluxOceanfloor size and grid
    if isnothing(rj.domain_oceanfloor.grid)
        @info "  creating oceanfloor grid, setting fluxOceanfloor size and grid"
        rj.domain_oceanfloor.grid = PB.Grids.UnstructuredVectorGrid(
            ncells=length(isedsurf)
        )
    else
        @info "  using existing oceanfloor.grid"
    end

    rj.domain_fluxOceanfloor  = PB.get_domain(model, "fluxOceanfloor")
    !isnothing(rj.domain_fluxOceanfloor) ||
        error("SedimentTransport.set_model_geometry no fluxOceanfloor Domain")
    if isnothing(rj.domain_fluxOceanfloor.grid)
        @info "  setting fluxOceanfloor grid"        
        rj.domain_fluxOceanfloor.grid = rj.domain_oceanfloor.grid
    else
        @info "  using existing fluxOceanfloor.grid"
    end

    PB.Grids.set_subdomain!(
        rj.domain_oceanfloor.grid,
        "sediment",
        PB.Grids.InteriorSubdomain(rj.domain.grid.ncells, isedsurf),
        true
    )

    # set sedimentfloor size and grid
    rj.domain_sedimentfloor = PB.get_domain(model, "sedimentfloor")
    !isnothing(rj.domain_sedimentfloor) ||
        error("SedimentTransport.set_model_geometry no sedimentfloor Domain")
    rj.domain_sedimentfloor.grid = PB.Grids.UnstructuredVectorGrid( ncells=length(isedfloor))
    PB.Grids.set_subdomain!(
        rj.domain_sedimentfloor.grid,
        "sediment",
        PB.Grids.InteriorSubdomain(rj.domain.grid.ncells, isedfloor),
        true
    )

    # if not already set by Ocean component, set fluxOceanBurial size and grid
    rj.domain_fluxOceanBurial = PB.get_domain(model, "fluxOceanBurial")
    !isnothing(rj.domain_fluxOceanBurial) ||
        error("SedimentTransport.set_model_geometry no fluxOceanBurial Domain")
    if isnothing(rj.domain_fluxOceanBurial.grid)
        @info "  setting fluxOceanBurial size and grid"
        rj.domain_fluxOceanBurial.grid = rj.domain_sedimentfloor.grid
    else
        @info "  using existing fluxOceanBurial.grid"
    end

    return nothing
end



function PB.register_methods!(rj::ReactionSedimentGridn1D)
    PB.add_method_setup!(
        rj,
        do_sediment_setup_grid,
        (PB.VarList_namedtuple(oceanfloor_grid_vars), PB.VarList_namedtuple(grid_vars)),
    )

    return nothing
end



function _calc_vertical_grid_linear(ncells, L)
    binedges = collect(range(0.0, -L, length=ncells+1))
    return (binedges[1:end-1], binedges[2:end], 0.5*(binedges[1:end-1]+binedges[2:end]))
end

"calculate a grid with quadratic spacing for z << eta, linear spacing for z >> eta"
function _calc_vertical_grid_quadratic(ncells, L, eta)

    # linear grid, NB: -L to get +ve values
    lin_zupper, lin_zlower, lin_zmid = _calc_vertical_grid_linear(ncells, -L)

    # transform z_lin (+ve) to z_xform (+ve)
    function z_xform(z_lin, L, eta)
        return L*(sqrt(z_lin^2 + eta^2) - eta)/(sqrt(L^2 + eta^2) - eta)
    end


    return (-z_xform.(lin_zupper, L, eta),
            -z_xform.(lin_zlower, L, eta),
            -z_xform.(lin_zmid, L, eta))
end

"set grid variables for 1D columns"
function do_sediment_setup_grid(
    m::PB.ReactionMethod,
    pars,
    (oceanfloor_grid_vars, grid_vars,),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    rj = m.reaction

    attribute_name == :setup || return

    if pars.f_grid[] == "linear"
        col_zupper, col_zlower, col_zmid = _calc_vertical_grid_linear(
            pars.ncellspercol[], pars.L[]
        )
    elseif pars.f_grid[] == "quadratic"
        col_zupper, col_zlower, col_zmid = _calc_vertical_grid_quadratic(
            pars.ncellspercol[], pars.L[], pars.grid_eta[]
        )
    else
        error("do_sediment_setup_grid $(PB.fullname(rj)): unknown f_grid=$(pars.f_grid[])")
    end


    for icol in 1:length(rj.domain.grid.Icolumns)
        colindices = rj.domain.grid.Icolumns[icol]  # indices for this column

        grid_vars.zupper[colindices]     .= col_zupper
        grid_vars.zlower[colindices]     .= col_zlower
        grid_vars.zmid[colindices]       .= col_zmid

        grid_vars.Abox[colindices]       .= oceanfloor_grid_vars.oceanfloor_Afloor[icol]

        grid_vars.pressure[colindices]   .= -1.0 .*(oceanfloor_grid_vars.oceanfloor_zfloor[icol] .+ grid_vars.zmid[colindices]) # dbar ~ depth in m

    end

    # attach coordinates to grid for output visualisation etc
    empty!(rj.domain.grid.z_coords)
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zmid", grid_vars.zmid, PB.get_variable(m, "zmid").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zlower", grid_vars.zlower, PB.get_variable(m, "zlower").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zupper", grid_vars.zupper, PB.get_variable(m, "zupper").attributes))
    
    grid_vars.volume .= grid_vars.Abox.*(grid_vars.zupper .- grid_vars.zlower)
    grid_vars.volume_total[] = sum(grid_vars.volume)

    grid_vars.rho_ref .= pars.rho_ref[] # kg m-3 assumed constant sw density conversion factor

    @info "do_sediment_setup_grid $(PB.fullname(rj)): volume_total $(grid_vars.volume_total[]) m^3"
    isfinite(grid_vars.volume_total[]) || error("configuration error - sediment volume_total is not finite")

    return nothing
end


end # module
