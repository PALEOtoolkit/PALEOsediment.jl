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

model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_FOAM_minimal_cfg.yaml"), 
    "sediment_Corg_Stest",
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        # oceanfloor oxygen
        # ("initial_value", "oceanfloor.O2_conc",  [0.200,    0.001,      0.200,      0.001] ),

        # bio rates
        # ("initial_value", "oceanfloor.alpha",  [465.0,     465.0,      114.0,      114.0]), # yr-1 bioirrigation coefficient at surface
        # ("initial_value", "oceanfloor.alpha",  0.5.*[465.0,     465.0,      114.0,      114.0]), # yr-1 reduced bioirrigation coefficient at surface
        # ("initial_value", "oceanfloor.alpha",  0.0), # yr-1 no bioirrigation 
        # ("initial_value", "oceanfloor.Dbio",  0.0), # yr-1 no bioturbation

    ]
) 

# tspan=(0.0, 1e5)
tspan=(0.0, 1e6)
# tspan=(0.0, 10.0)

initial_state, modeldata = PALEOmodel.initialize!(model)

run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB requires min/max ratio for robustness so check all tracers initial_value is not zero
# newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, 1e6, 0.1, 10.0
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-30, 1e6, 1e-2, 1e2
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    run, initial_state, modeldata, tspan, 1e-6,
    deltat_fac=2.0,
    solvekwargs=(
        ftol=1e-7,
        # ftol=1e-3,
        # iterations=20,
        iterations=50, # CFA formation has a discontinuous derivative -> need more iterations ?
        method=:newton,
        linesearch=LineSearches.Static(),
        apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(newton_min, newton_max, newton_min_ratio, newton_max_ratio),
        linsolve=PALEOmodel.SolverFunctions.SparseLinsolveUMFPACK(),
        # store_trace=true,
        # show_trace=true, 
    ),
    verbose=false,
) 


############################
# Plot 
############################

# single plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# multiple plots per screen
colrange=1:4
gr(size=(1200, 800))
pager = PALEOmodel.PlotPager((1, 4), (legend_background_color=nothing, margin=(5, :mm)))

solute_burial_flux=true
plot_phi(run.output; colrange, pager)
plot_w(run.output; colrange, pager)
plot_biorates(run.output; colrange, pager)
plot_Corg_O2(run.output; Corgs=["Corg",], colT=[first(tspan), last(tspan)], colrange, pager)
plot_Corg_RCmultiG(run.output; Corg_indices=1:12, colrange, pager)
plot_solutes(run.output; colT=[first(tspan), last(tspan)], solutes=["P", "DIC", "TAlk",  "SO4", "H2S", "CH4", ], colrange, pager)
plot_budget(run.output; name="C", solids=["Corg", ], solutes=["DIC", "CH4"], solute_burial_flux, pager, colrange, fluxylims=(-4, 4), concxscale=:log10, concxlims=(1e-3, 1e4))
plot_budget(run.output; name="P", solids=["Corg"], stoich_factors=Dict("Corg"=>1/106), solutes=["P"], solute_burial_flux, pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e3))
plot_budget(run.output; name="S", solids=[], solutes=["H2S", "SO4"], solute_burial_flux, pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4), fluxylims=(-1, 1) )
plot_rates(run.output; colT=[first(tspan), last(tspan)], plot_freminOrgTot=false, remin_rates=["reminOrgOxO2", "reminOrgOxSO4", "reminOrgOxCH4"], colrange, pager)
plot_conc_summary(run.output; species=["DIC", "O2", "P", "SO4", "H2S", "CH4"], pager, colrange, xscale=:log10, xlims=(1e-3, 1e2))
pager = PALEOmodel.PlotPager((2, 4), (legend_background_color=nothing, margin=(5, :mm)))
plot_summary_stacked(run.output; species=["SO4", "CH4", "DIC", "P", "O2", "H2S"], plot_pHfree=false, pager, colrange)
pager(:newpage)


