using Logging
using Sundials
import LineSearches
using BenchmarkTools

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOaqchem
import PALEOsediment

###############################################################
# 10 sediment columns with oceanfloor [O2] and [SO4] gradients
########################################################

global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_sediment_expts.jl")
include("../plot_sediment.jl")


model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_cfg.yaml"), 
    "sediment_Corg_O2";
    modelpars=Dict("num_columns"=>10),
)

#############################
# Steady state solutions
############################
    
config_sediment_expts(model, 
    [
        "baseline",
        # 2 component Corg input
        ("initial_value", "fluxOceanfloor.particulateflux_Corg1", 130e-2), # mol yr-1 = umol/cm^2/y * 1e-6 * area m^2 * 1e4
        ("initial_value", "fluxOceanfloor.particulateflux_Corg2", 18.5e-2), # mol yr-1 = umol/cm^2/y * 1e-6 * area m^2 * 1e4
        # physical pars for 10 "shelf" columns
        ("initial_value", "oceanfloor.zfloor", -100.0), # m
        ("initial_value", "oceanfloor.temp", 278.15), # K
        ("initial_value", "oceanfloor.w_accum", 0.03e-2), # m yr-1
        ("initial_value", "oceanfloor.Dbio", 0.0), # m^2 yr-1
        # oceanfloor solute boundary conditions
        ("initial_value", "oceanfloor.O2_conc", [0.250,  0.125, 0.050, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), # mol m-3
        ("initial_value", "oceanfloor.SO4_conc",  [28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, 10000e-3, 5000e-3, 1000e-3, 100e-3, 0.0]), # mol m-3
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
gr(size=(1500, 900))
pager = PALEOmodel.PlotPager((2,5), (legend_background_color=nothing, margin=(5, :mm)))

colrange=1:10 # number of columns
plot_Corg_O2(run.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], colrange, pager=pager)
plot_solutes(run.output; colT=[first(tspan), last(tspan)], colrange, pager=pager)
plot_rates(run.output; colT=[first(tspan), last(tspan)], colrange, pager=pager)

pager(:newpage)


# summary plots of solute fluxes vs O2_conc
O2_conc = [PALEOmodel.get_array(run.output, "oceanfloor.O2_conc", (tmodel=1e12, cell=i)).values for i in 1:10]
SO4_conc = [PALEOmodel.get_array(run.output, "oceanfloor.SO4_conc", (tmodel=1e12, cell=i)).values for i in 1:10]
soluteflux_O2 = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_O2", (tmodel=1e12, cell=i)).values for i in 1:10]
soluteflux_H2S = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_H2S", (tmodel=1e12, cell=i)).values for i in 1:10]
soluteflux_CH4 = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_CH4", (tmodel=1e12, cell=i)).values for i in 1:10]

reminOrgOxO2_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxO2", (tmodel=1e12, column=i)).values) for i in 1:10]
reminOrgOxSO4_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxSO4", (tmodel=1e12, column=i)).values) for i in 1:10]
reminOrgOxCH4_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxCH4", (tmodel=1e12, column=i)).values) for i in 1:10]

redox_H2S_O2_total = [sum(PALEOmodel.get_array(run.output, "sediment.redox_H2S_O2", (tmodel=1e12, column=i)).values) for i in 1:10]

gr(size=(800, 600))
pager = PALEOmodel.PlotPager((2,2), (legend_background_color=nothing, margin=(5, :mm)))
pager(
    (
        p = plot(title="solute fluxes [SO4] = 28 mM", xlabel="[O2] (mol m-3)", ylabel="flux (mol m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true);
        plot!(p, O2_conc[1:5], soluteflux_O2[1:5]; label="O2");
        plot!(p, O2_conc[1:5], soluteflux_H2S[1:5]; label="H2S");
        plot!(p, O2_conc[1:5], soluteflux_CH4[1:5]; label="CH4");
    ),
    (
        p = plot(title="solute fluxes [O2] = 0", xlabel="[SO4] (mol m-3)", ylabel="flux (mol m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true);
        plot!(p, SO4_conc[6:10], soluteflux_O2[6:10]; label="O2");
        plot!(p, SO4_conc[6:10], soluteflux_H2S[6:10]; label="H2S");
        plot!(p, SO4_conc[6:10], soluteflux_CH4[6:10]; label="CH4");
    ),
    (
        p = plot(title="remin fluxes [SO4] = 28 mM", xlabel="[O2] (mol m-3)", ylabel="flux (mol O2eq m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true);
        plot!(p, O2_conc[1:5], reminOrgOxO2_total[1:5]; label="O2");
        plot!(p, O2_conc[1:5], reminOrgOxSO4_total[1:5]; label="SO4");
        plot!(p, O2_conc[1:5], reminOrgOxCH4_total[1:5]; label="CH4");
        plot!(p, O2_conc[1:5], redox_H2S_O2_total[1:5]; label="redox H2S + O2");
    ),
   
    (
        p = plot(title="remin fluxes [O2] = 0", xlabel="[SO4] (mol m-3)", ylabel="flux (mol O2eq m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true);
        plot!(p, SO4_conc[6:10], reminOrgOxO2_total[6:10]; label="O2");
        plot!(p, SO4_conc[6:10], reminOrgOxSO4_total[6:10]; label="SO4");
        plot!(p, SO4_conc[6:10], reminOrgOxCH4_total[6:10]; label="CH4");        
    ),
    :newpage,
)