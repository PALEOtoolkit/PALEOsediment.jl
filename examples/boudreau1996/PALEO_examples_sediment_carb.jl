using Logging
using Sundials
import LineSearches
using BenchmarkTools

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOaqchem
import PALEOsediment

global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_sediment_expts.jl")
include("../plot_sediment.jl")


model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_cfg.yaml"), 
    "sediment_Corg_O2_carb",
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        "baseline",
    ]
) 

tspan=(0.0, 10000.0)

initial_state, modeldata = PALEOmodel.initialize!(model)

run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())


# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB some tracers start at zero so set newton_max_ratio=Inf
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, Inf, 0.1, Inf 
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    run, initial_state, modeldata, tspan, 1e-3,
    deltat_fac=2.0,
    solvekwargs=(
        ftol=1e-7,
        iterations=20,
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
pager = PALEOmodel.PlotPager((1,3), (legend_background_color=nothing, ))

plot_Corg_O2(run.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], pager)
plot_solutes(run.output; solutes = ["P", "SO4", "H2S", "CH4", "DIC", "TAlk"], colT=[first(tspan), last(tspan)], pager)
plot_rates(run.output; colT=[first(tspan), last(tspan)], pager=pager)
plot_carbchem(run.output; include_constraint_error=true, colT=last(tspan), pager)

pager(:newpage)


