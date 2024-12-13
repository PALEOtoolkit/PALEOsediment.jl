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
    joinpath(@__DIR__, "PALEO_transport_RCmultiG_cfg.yaml"), 
    "sediment_abiotic_O2",
    modelpars=Dict("CIsotope"=>"IsotopeLinear"),  # enable C isotopes
)

#############################
# Steady state solutions
############################

tspan=(0.0, 1e5)
toutput=[0.0, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5]

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB some tracers start at zero so set newton_max_ratio=Inf
# TODO isotope mol*delta can be -ve so can't use StepClampMultAll! for Newton robustness 
# newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, Inf, 0.1, Inf 
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-3;
    deltat_fac=2.0,
    saveat=toutput,
    solvekwargs=(
        # ftol=1e-7,
        ftol=1e-6,
        iterations=20,
        method=:newton,
        linesearch=LineSearches.Static(),
        # apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(newton_min, newton_max, newton_min_ratio, newton_max_ratio),
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
gr(size=(1200, 800))
pager = PALEOmodel.PlotPager((1,4), (legend_background_color=nothing, margin=(5, :mm)))
colrange=1:4 # number of columns

plot_phi(paleorun.output; colrange, pager)
plot_w(paleorun.output; colrange, pager)
plot_biorates(paleorun.output; colrange, pager)
plot_Corg_O2(paleorun.output, Corgs=["Corg"]; colT=toutput, colrange, pager)
plot_solutes(paleorun.output; solutes=["DIC"], colT=toutput, colrange, pager)
plot_solid_deltas(paleorun.output; solids=["Corg"], colT=toutput, colrange, pager)
plot_solute_deltas(paleorun.output; solutes=["DIC"], colT=toutput, colrange, pager)
plot_Corg_RCmultiG(paleorun.output; colrange, pager)

pager(:newpage)


