using Logging
using Sundials
import LineSearches
using BenchmarkTools

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOsediment


global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_sediment_expts.jl")
include("../plot_sediment.jl")
include("../SumColumns_dev.jl")

model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_FOAM_cfg.yaml"), 
    "sediment_FOAM",
)

# set true for additional plots
verbose_plots = false

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        # oceanfloor oxygen
        # ("initial_value", "oceanfloor.O2_conc",  [0.200,    0.001,      0.200,      0.001] ),

        # bio rates
        # ("initial_value", "oceanfloor.alpha",  0.0), # yr-1 no bioirrigation 
        # ("initial_value", "oceanfloor.Dbio",  0.0), # yr-1 no bioturbation

        # FeII adsorption
        # Zhao eqns Table S4  quote K ~ Fe adsorbed / Fe solute ~ 500 !! (without corrections for porosity of ~2x)
        # PALEO defines  [FeIItot] * volume = [FeIIsolute]*phi*volume + K_eqb * [FeIIsolute] *(1-phi)*volume
        #              so Fe in solute / Fe adsorbed = phi / (K_eqb  * (1 - phi))
        #                                            = 1/K_eqb * phi / (1 - phi)
        #                  Fe adsorbed / Fe in solute = K_eqb * (1 - phi) / phi
        #
        # ("set_par", "sediment", "FeIIadsorb", "K_eqb", 500.0), # Zhao (2020) default ?! (or units problem?)
        ("set_par", "sediment", "FeIIadsorb", "K_eqb", 15.0), # 
        # ("set_par", "sediment", "FeIIadsorb", "K_eqb", 0.0),  # no adsorption

        # biotite input: 0.04 mol m-2 yr-1 = 0.004 mmol cm-2 yr-1 is the Zhao 2020 value
        # NB: no actual chemical biotite decay, constant biotite with depth set by input vs transport
        ("initial_value", "fluxOceanfloor.particulateflux_Biotite", [0.025, 0.03, 0.035, 0.04]),

        # reduce aragonite & calcite precipitation rate
        # ("set_par", "sediment", "CaCO3arag_precipdissol", "K_precip", 0.75e-3), # 3e-3 yr-1 Zhao (2020)
        # ("set_par", "sediment", "CaCO3calc_precipdissol", "K_precip", 0.75e-4), # 3e-4 yr-1 Zhao (2020)
        ("set_par", "sediment", "CaCO3arag_precipdissol", "K_precip", 1.5e-3), # 3e-3 yr-1 Zhao (2020)
        ("set_par", "sediment", "CaCO3calc_precipdissol", "K_precip", 1.5e-4), # 3e-4 yr-1 Zhao (2020)


    ]
) 

tspan=(0.0, 1e5)
# tspan=(0.0, 1e6)
# tspan=(0.0, 1e0)
# tspan=(0.0, 1e2)

initial_state, modeldata = PALEOmodel.initialize!(model)

run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB requires min/max ratio for robustness so check all tracers initial_value is not zero
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, 1e6, 0.1, 10.0
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    run, initial_state, modeldata, tspan, 1e-6,
    deltat_fac=2.0,
    solvekwargs=(
        # ftol=1e-7,
        ftol=1e-3,
        iterations=50, 
        method=:newton,
        linesearch=LineSearches.Static(),
        apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(newton_min, newton_max, newton_min_ratio, newton_max_ratio),
        # store_trace=true,
        # show_trace=true, 
    ),
    verbose=false,
) 


############################
# Plot 
############################

# multiple plots per screen
colrange=1:4
gr(size=(1200, 800)) # screen size

solute_burial_flux=true
if verbose_plots    
    pager = PALEOmodel.PlotPager((1, 4), (legend_background_color=nothing, margin=(5, :mm)))

    plot_phi(run.output; colrange, pager)
    plot_w(run.output; colrange, pager)
    plot_biorates(run.output; colrange, pager)
    plot_Corg_O2(run.output; Corgs=["Corg",], colT=[first(tspan), last(tspan)], colrange, pager)
    plot_Corg_RCmultiG(run.output; Corg_indices=1:12, colrange, pager)
    plot_solutes(run.output; colT=[first(tspan), last(tspan)], solutes=["TP", "DIC", "TAlk",  "NO3", "TNH3", "SO4", "TH2S", "CH4",  "H2", "MnII",], colrange, pager)
    # plot_sediment_FeS_summary(run.output; colrange, pager)
    plot_solids(run.output; colT=[first(tspan), last(tspan)], solids=["MnHR", "MnMR", "FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr",], colrange, pager)
    plot_carbchem(run.output; include_constraint_error=true, colT=last(tspan), colrange, pager)
    plot_budget(run.output; name="C", solids=["Corg", "CaCO3calc", "CaCO3arag", "MnCO3rhod"], solutes=["DIC", "CH4"], solute_burial_flux,
        pager, colrange, fluxylims=(-4, 4), concxscale=:log10, concxlims=(1e-3, 1e4))
    plot_budget(run.output; name="TP", solids=["Corg", "PFeHR", "PFeMR", "CFA"], stoich_factors=Dict("Corg"=>1/106), solutes=["TP"], solute_burial_flux,
        pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e3))
    plot_budget(run.output; name="Mn", solids=["MnHR", "MnMR", "MnCO3rhod"], solutes=["MnII"], solute_burial_flux,
        pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4))
    # don't include solute_burial_flux as FeIItot is already counted once as a solid
    plot_budget(run.output; name="Fe", solids=["FeHR", "FeMR", "FePR", "FeMag", "FeSm", "FeS2pyr", "FeIItot"], solutes=["FeIItot"], extras=["budgets.Biotite_dissolflux"], solute_burial_flux=false,
        stoich_factors=Dict("FeMag"=>3, "budgets.Biotite_dissolflux"=>2),
        pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4))
    plot_budget(run.output; name="S", solids=["FeSm", "FeS2pyr", "S0"], solutes=["TH2S", "SO4"], solute_burial_flux,
        stoich_factors=Dict("FeS2pyr"=>2,),
        pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4), fluxylims=(-1, 1) )
    # TODO TAlk conservation requires counting adsorbed FeII that is buried ?!
    plot_budget(run.output; name="TAlk", solids=["FeIItot"], solutes=["TAlk", "SO4", "NO3", "TNH3", "TP", "MnII", "FeIItot", "Ca", "Mg", "K", "Na", "F"], solute_burial_flux,
        stoich_factors=Dict("SO4"=>2, "NO3"=>1, "TNH3"=>-1, "TP"=>1, "MnII"=>-2, "FeIItot"=>-2, "Ca"=>-2, "Mg"=>-2, "K"=>-1), 
        pager, colrange, concxscale=:identity, concxlims=(-Inf, Inf), fluxylims=(-10, 10) )
    plot_rates(run.output; colT=[first(tspan), last(tspan)], plot_freminOrgTot=false, remin_rates=["reminOrgOxO2", "reminOrgOxNO3", "reminOrgOxSO4", "reminOrgOxCH4"], colrange, pager)
    plot_conc_summary(run.output; species=["DIC", "O2", "TP", "SO4", "TH2S", "CH4"], pager, colrange, xscale=:log10, xlims=(1e-3, 1e2))
end

# summary plots
pager = PALEOmodel.PlotPager((2, 4), (legend_background_color=nothing, margin=(5, :mm)))
plot_budgets(
    run.output;
    budgets="budgets.net_input_".*["C", "P", "S", "Mn", "Fe", "TAlk"],
    ylims=(-1e-6, 1e-6),
    pager, colrange,
)
pager(:newpage)

plot_summary_stacked(
    run.output; 
    species=["MnII", "FeII", "SO4", "CH4", "FeMag", "TNH3", "Ca", "DIC", "TP", "F", "CFA", "O2", "NO3", "TH2S", "TAlk",],
    others=["NotOmega_CFA",],
    pager, colrange
)
pager(:newpage)


