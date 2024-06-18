using Logging
using Sundials
import LineSearches
using BenchmarkTools

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOsediment

do_benchmarks = false

global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_sediment_expts.jl")
include("../plot_sediment.jl")
include("../SumColumns_dev.jl")

model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_ironsulphur_cfg.yaml"), 
    "sediment_Corg_O2_Fe",
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        "baseline",
    ]
) 

# tspan=(0.0, 10000.0)
tspan=(0.0, 1e6)

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB some tracers start at zero so set newton_max_ratio=Inf
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, Inf, 0.1, Inf 
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-3,
    deltat_fac=2.0,
    solvekwargs=(
        ftol=1e-7,
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

# single plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# multiple plots per screen
gr(size=(900, 600))
pager = PALEOmodel.PlotPager((1,3), (legend_background_color=nothing, margin=(5, :mm)))
colrange=1:3

plot_budgets(
    paleorun.output;
    budgets="budgets.net_input_".*["C", "P", "S", "Fe", "TAlk"],
    ylims=(-1e-6, 1e-6),
    pager, colrange,
)
pager(:newpage)

plot_Corg_O2(paleorun.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], colrange, pager)
plot_solutes(paleorun.output; colT=[first(tspan), last(tspan)], solutes=["P", "SO4", "H2S", "CH4", "FeII"], colrange, pager)
plot_solids(paleorun.output; colT=[first(tspan), last(tspan)], solids=["FeHR", "FeMR", "FePR"], colrange, pager)
plot_rates(paleorun.output; colT=[first(tspan), last(tspan)], remin_rates=["reminOrgOxO2", "reminOrgOxFeIIIOx", "reminOrgOxSO4", "reminOrgOxCH4"], colrange, pager)

pager(:newpage)


