using Logging
using Sundials
import LineSearches
using BenchmarkTools

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOaqchem
import PALEOsediment

do_benchmarks = false

global_logger(ConsoleLogger(stderr,Logging.Info))

include("../plot_sediment.jl")


model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_transport_mud_cfg.yaml"), 
    "sediment_abiotic_O2",
)

#############################
# Steady state solutions
############################

tspan=(0.0, 1e4)
toutput=[0.0, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB some tracers start at zero so set newton_max_ratio=Inf
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, Inf, 0.1, Inf 
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-3;
    deltat_fac=2.0,
    saveat=toutput,
    solvekwargs=(
        ftol=1e-7,
        iterations=20,
        method=:newton,
        linesearch=LineSearches.Static(),
        apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(newton_min, newton_max, newton_min_ratio, newton_max_ratio),
        always_step=true,
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

plot_w(paleorun.output; pager=pager)
plot_solids(paleorun.output, solids=["M"]; colT=toutput, pager=pager)
plot_solids_volume_frac(paleorun.output, solids=["M"]; colT=toutput, pager=pager)
plot_solutes(paleorun.output, solutes=["O2"]; colT=toutput, pager=pager)

pager(:newpage)


