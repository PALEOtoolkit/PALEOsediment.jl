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

include("config_sediment_expts.jl")
include("../plot_sediment.jl")


model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_cfg.yaml"), 
    "sediment_Corg_O2",
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

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# Newton, no line search
# sol = PALEOmodel.SteadyState.steadystateForwardDiff(
#     paleorun, initial_state, modeldata, 0.0,
#     solvekwargs=(
#         ftol=5e-9,
#         method=:newton,
#         store_trace=true,
#         show_trace=true
#     )
# );

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB some tracers start at zero so set newton_max_ratio=Inf
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, Inf, 0.1, Inf 
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-3,
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

##########################
# Time dependent solutions
##########################
#=
config_sediment_expts(model, 
    [
        "baseline",
        "slowredox",
    ]
) 
tspan=(0,10000.0)
initial_state, modeldata = PALEOmodel.initialize!(model)
paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())
sol = PALEOmodel.ODE.integrateForwardDiff(paleorun, initial_state, modeldata, tspan, 
#    solvekwargs=(reltol=1e-3, abstol=1e-4*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
#    solvekwargs=(reltol=1e-3, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
    solvekwargs=(reltol=1e-5, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
#    solvekwargs=(reltol=1e-7, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
#    solvekwargs=(reltol=1e-9, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
#    solvekwargs=(reltol=1e-11, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=10000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
#    solvekwargs=(reltol=1e-13, abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), maxiters=1000000, saveat=[0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]))  # first run includes JIT time
# sol = PALEOmodel.ODE.integrateDAEForwardDiff(paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-5,))  # first run includes JIT time

# Test additional solver options
# sol = PALEOmodel.ODE.integrateDAEForwardDiff(paleorun, initial_state, modeldata, tspan, jac_ad=:ForwardDiff, alg=IDA(), solvekwargs=(reltol=1e-5,))
# sol = PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-5,))  # first run includes JIT time
=#


############################
# Plot 
############################

# single plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# multiple plots per screen
gr(size=(900, 600))
pager = PALEOmodel.PlotPager((1,3), (legend_background_color=nothing, margin=(5, :mm)))

plot_Corg_O2(paleorun.output, Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], pager=pager)
plot_solutes(paleorun.output, colT=[first(tspan), last(tspan)], pager=pager)
plot_rates(paleorun.output, colT=[first(tspan), last(tspan)], pager=pager)

pager(:newpage)


