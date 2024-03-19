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
        # ("set_par", "sediment", "PyrH2S", "R_Pyr_H2S", 0.0), # zero rate (disable) pyrite formation
        ("set_par", "sediment", "PyrH2S", "R_Pyr_H2S", 1e2), # 1e5 M yr-1 = 1e5*1e-3 (mol m-3) yr-1, Dale (2015)
        # Pyrite oxidation rate
        ("set_par", "sediment", "redox_FeS2pyr_O2", "R_FeS2pyr_O2",  1.0),  # (mol m-3)-1 yr-1,  1e3 M-1 yr-1 = 1e3*1e-3, Dale (2015)
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "R_FeS2pyr_O2",  0.0), # disable pyrite oxidation
        #
        # no Fe input
        # ("initial_value", "fluxOceanfloor.particulateflux_FeHR", 0.0), 
        # ("initial_value", "fluxOceanfloor.particulateflux_FeMR", 0.0), 
        # ("initial_value", "fluxOceanfloor.particulateflux_FePR", 0.0), 
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
pager = PALEOmodel.PlotPager((1,3), (legend_background_color=nothing, margin=(5, :mm)))

plot_Corg_O2(run.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], pager=pager)
plot_solutes(run.output; colT=[first(tspan), last(tspan)], solutes=["P", "SO4", "SmIIaqtot", "CH4", "H2", "FeIIaqtot", "DIC", "TAlk"], pager=pager)
plot_sediment_FeS_summary(run.output; FeII_species = ["FeIIaqtot", "FeII", "FeSaq",], pager=pager)
plot_solids(run.output; colT=[first(tspan), last(tspan)], solids=["FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr"], pager=pager)
plot_rates(run.output; colT=[first(tspan), last(tspan)], remin_rates=["reminOrgOxO2", "reminOrgOxFeIIIOx", "reminOrgOxSO4", "reminOrgOxCH4"], pager=pager)
plot_carbchem(run.output; include_constraint_error=true, colT=last(tspan), pager)

pager(:newpage)


