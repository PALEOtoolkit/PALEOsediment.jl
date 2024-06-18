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
include("CFA.jl")

model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_NNFeMnP_cfg.yaml"), 
    "sediment_Corg_O2NNMnFeSP",
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        # oceanfloor oxygen
        ("initial_value", "oceanfloor.O2_conc",  [0.200,    0.001,      0.200,      0.001] ),

        # Corg input
        # sensitivity test: remove most-reactive bin, maintaining same Corg input
        # ("set_par", "sediment", "reservoir_Corg", "k_dist_modifier", [1.0,   1.0, 1.0,   1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, 1.0,  1.0, 0.0]), # remove 1 highest-k fractions

        # Fe & Mn input and remin
        # high value 75 umol cm-2 yr-1 (0.75 mol m-2 yr-1, Van Cappellen 1996 test case)
        # ("initial_value", "fluxOceanfloor.particulateflux_FeHR", 7.0.*[0.112, 0.112,  0.0177, 0.0177]), # mol m-2 yr-1
        # high value 40 um cm-2 yr-1 (0.4 mol m-2 yr-1, Van Cappellen 1996 test case)
        # ("initial_value", "fluxOceanfloor.particulateflux_MnHR", 12.0.*[0.0342, 0.0342, 0.0101, 0.0101]), # mol m-2 yr-1
        # low values for limiting conc
        # ("set_par", "sediment", "reminsed", "MnIVOxreminlimit", 40.0), # (mol m-3) default = 16e-6*2.5e6 (16 umol MnO2 g-1 assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3), Van Cappellen 1996)
        # ("set_par", "sediment", "reminsed", "FeIIIOxreminlimit", 250.0), # (mol m-3) default = 100e-6*2.5e6 (100 umol Fe(OH)3 g-1 assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3), Van Cappellen 1996)
        # ("set_par", "sediment", "reminsed", "FeIIIOxreminlimit", 0.1*250.0), # (mol m-3) = 100e-6*2.5e6 (100 umol Fe(OH)3 g-1 assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3), Van Cappellen 1996)


        # bio rates
        ("initial_value", "oceanfloor.alpha",  [465.0,     465.0,      114.0,      114.0]), # yr-1 bioirrigation coefficient at surface
        # ("initial_value", "oceanfloor.alpha",  0.5.*[465.0,     465.0,      114.0,      114.0]), # yr-1 reduced bioirrigation coefficient at surface
        # ("initial_value", "oceanfloor.alpha",  0.0), # yr-1 no bioirrigation 

        # Secondary redox rates         
        # ("set_par", "sediment", "redox_H2S_O2", "K",  1e2),  # (mol m-3)-1 yr-1,  1e5 M-1 yr-1 -> 1e3*1e-3, Dale (2015) R15
        # ("set_par", "sediment", "redox_H2S_O2", "K",  1e5),  # (mol m-3)-1 yr-1 sensitivity test to x 1e3 higher rate

        # Fe-S system
        # Pyrite formation rate
        # ("set_par", "sediment", "pyrite_H2S", "K", 0.0), # zero rate (disable) pyrite formation
        ("set_par", "sediment", "pyrite_H2S", "K", 1e2), # 1e5 M yr-1 = 1e5*1e-3 (mol m-3) yr-1, Dale (2015)
        # pyrite oxidation rate
        ("set_par", "sediment", "redox_FeS2pyr_O2", "K",  1.0),  # (mol m-3)-1 yr-1,  1e3 M-1 yr-1 = 1e3*1e-3, Dale (2015)
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "K",  100.0),  # test x100 rate
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "K",  0.0), # disable pyrite oxidation
    ]
) 

tspan=(0.0, 100000.0)
# tspan=(0.0, 10.0)

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# PTC, Newton, no line search
# Bounds and max step size for Newton solve. NB requires min/max ratio for robustness so check all tracers initial_value is not zero
# newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-80, 1e6, 0.1, 10.0
newton_min, newton_max, newton_min_ratio, newton_max_ratio = 1e-30, 1e6, 1e-2, 1e2
PALEOmodel.SteadyState.steadystate_ptcForwardDiff(
    paleorun, initial_state, modeldata, tspan, 1e-6,
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

plot_phi(paleorun.output; colrange, pager)
plot_w(paleorun.output; colrange, pager)
plot_biorates(paleorun.output; colrange, pager)
plot_Corg_O2(paleorun.output; Corgs=["Corg",], colT=[first(tspan), last(tspan)], colrange, pager)
plot_Corg_RCmultiG(paleorun.output; colrange, pager)
plot_solutes(paleorun.output; colT=[first(tspan), last(tspan)], solutes=["TP", "DIC", "NO3", "NO2", "TNH3", "SO4", "TH2S", "CH4", "H2", "MnII", "FeII"], colrange, pager)
plot_solids(paleorun.output; colT=[first(tspan), last(tspan)], solids=["MnHR", "MnMR", "FeHR", "FeMR", "FePR", "S0", "FeSm", "FeS2pyr",], colrange, pager)
plot_budget(paleorun.output; name="P", solids=["Corg", "PFeHR", "PFeMR", "CFA"], stoich_factors=Dict("Corg"=>1/106), solutes=["TP"], solute_burial_flux,
    pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e3))
plot_budget(paleorun.output; name="Mn", solids=["MnHR", "MnMR"], solutes=["MnII"], solute_burial_flux,
    pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e3))
plot_budget(paleorun.output; name="Fe", solids=["FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr"], solutes=["FeII"], solute_burial_flux,
    pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4))
plot_budget(paleorun.output; name="S", solids=["S0", "FeSm", "FeS2pyr"], solutes=["TH2S", "SO4"], stoich_factors=Dict("FeS2pyr"=>2), solute_burial_flux,
    pager, colrange, concxscale=:log10, concxlims=(1e-3, 1e4), fluxylims=(-1, 1) )
plot_rates(paleorun.output; colT=[first(tspan), last(tspan)], remin_rates=["reminOrgOxO2", "reminOrgOxNO2", "reminOrgOxNO3NO2", "reminOrgOxMnIVOx", "reminOrgOxFeIIIOx", "reminOrgOxSO4", "reminOrgOxCH4"], colrange, pager)

plot_budgets(
    paleorun.output;
    budgets="budgets.net_input_".*["C", "N", "P", "S", "Mn", "Fe", "TAlk"],
    ylims=(-1e-6, 1e-6),
    pager, colrange,
)
pager(:newpage)

plot_conc_summary(paleorun.output; species=["O2", "NO3", "NO2", "TNH3", "TP", "TH2S", "FeII", "MnII", "CH4", "F"], pager, colrange, xscale=:log10, xlims=(1e-3, 1e0))
plot_conc_summary(paleorun.output; species=["Corg", "FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr", "PFeHR", "PFeMR", "CFA"], pager, colrange, xscale=:log10, xlims=(1e-3, 1e4))

pager(:newpage)


