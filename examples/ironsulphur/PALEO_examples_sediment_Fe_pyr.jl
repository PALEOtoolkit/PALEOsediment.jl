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
    joinpath(@__DIR__, "PALEO_examples_sediment_ironsulphur_cfg.yaml"), 
    "sediment_Corg_O2_Fe_pyr_carb",
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        "baseline",
        # Fe-S system
        # Pyrite formation rate
        # ("set_par", "sediment", "pyrite_H2S", "K", 0.0), # zero rate (disable) pyrite formation
        ("set_par", "sediment", "pyrite_H2S", "K", 1e2), # 1e5 M yr-1 = 1e5*1e-3 (mol m-3) yr-1, Dale (2015)
        # Pyrite oxidation rate
        ("set_par", "sediment", "redox_FeS2pyr_O2", "K",  1.0),  # (mol m-3)-1 yr-1,  1e3 M-1 yr-1 = 1e3*1e-3, Dale (2015)
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "K",  0.0), # disable pyrite oxidation
        #
        # no Fe input
        # ("initial_value", "fluxOceanfloor.particulateflux_FeHR", 0.0), 
        # ("initial_value", "fluxOceanfloor.particulateflux_FeMR", 0.0), 
        # ("initial_value", "fluxOceanfloor.particulateflux_FePR", 0.0), 
    ]
) 

tspan=(0.0, 100000.0)

# tspan=(0.0, 100.0)

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB requires min/max ratio for robustness so check all tracers initial_value is not zero
# newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, 1e6, 0.1, 10.0
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-30, 1e6, 1e-2, 1e2
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-3,
    deltat_fac=2.0,
    solvekwargs=(
        ftol=1e-7,
        # ftol=1e-3,
        iterations=30,
        method=:newton,
        linesearch=LineSearches.Static(),
        apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(newton_min, newton_max, newton_min_ratio, newton_max_ratio),
        # linsolve=PALEOmodel.SolverFunctions.SparseLinsolveUMFPACK(),
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

plot_Corg_O2(paleorun.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], pager=pager)
plot_solutes(paleorun.output; colT=[first(tspan), last(tspan)], solutes=["P", "SO4", "TH2S", "CH4", "H2", "TFeII", "DIC", "TAlk"], pager=pager)
plot_sediment_FeS_summary(paleorun.output; FeII_species = ["TFeII", "FeII", "FeSaq",], pager=pager)
plot_solids(paleorun.output; colT=[first(tspan), last(tspan)], solids=["FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr"], pager=pager)
plot_rates(paleorun.output; colT=[first(tspan), last(tspan)], remin_rates=["reminOrgOxO2", "reminOrgOxFeIIIOx", "reminOrgOxSO4", "reminOrgOxCH4"], pager=pager)
plot_carbchem(paleorun.output; include_constraint_error=true, colT=last(tspan), pager)

pager(:newpage)


