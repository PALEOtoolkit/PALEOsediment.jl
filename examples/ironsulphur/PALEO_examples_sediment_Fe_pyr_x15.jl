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

###############################################################################
# Set of 15 columns, all with same shelf environment, Corg input with two reactivity fractions 
# (a combination of the high-reactivity fraction from the shelf/slope case, and the low reactivity
# fraction from the rise case), excess Fe input:
# - Columns 1-5 with [SO4] 28 mM (high modern value), decreasing values of [O2], with bioturbation
# - Columns 6-10 with [SO4] 28 mM (high modern value), decreasing values of [O2], no bioturbation
# - Columns 11-15 with zero [O2], no bioturbation decreasing values of [SO4]
##################################################################################

num_columns=15
model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_sediment_ironsulphur_cfg.yaml"), 
    "sediment_Corg_O2_Fe_pyr_carb";
    modelpars=Dict("num_columns"=>num_columns),
)

#############################
# Steady state solutions
############################

# Fe oxide input 1000 FeT umol/m^2/d, 1/6 to each of FeHR, FeMR, FePR
# FeT = 0.36525 # mol yr-1 = umol/m^2/d * 1e-6 * area m^2 * 365.25
# high value 75 umol cm-2 yr-1 FeHR (0.75 mol m-2 yr-1, Van Cappellen 1996 test case)
FeT = 10*0.36525

config_sediment_expts(model, 
    [
        # physical pars for 10 "shelf" columns
        ("initial_value", "oceanfloor.zfloor", -100.0), # m
        ("initial_value", "oceanfloor.temp", 278.15), # K
        ("initial_value", "oceanfloor.w_accum", 0.03e-2), # m yr-1
        ("initial_value", "oceanfloor.Dbio", 
            [
                1.35e-4, 1.35e-4, 1.35e-4, 1.35e-4, 0.0,  # 5 oxic bioturbated
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
            ]
        ), # m^2 yr-1

        # 2 component Corg input
        ("initial_value", "fluxOceanfloor.particulateflux_Corg1", 130e-2), # mol yr-1 = umol/cm^2/y * 1e-6 * area m^2 * 1e4
        ("initial_value", "fluxOceanfloor.particulateflux_Corg2", 18.5e-2), # mol yr-1 = umol/cm^2/y * 1e-6 * area m^2 * 1e4
        # Fe oxide input 1000 FeT umol/m^2/d, 1/6 to each of FeHR, FeMR, FePR
        ("initial_value", "fluxOceanfloor.particulateflux_FeHR", FeT/6), 
        ("initial_value", "fluxOceanfloor.particulateflux_FeMR", FeT/6), 
        ("initial_value", "fluxOceanfloor.particulateflux_FePR", FeT/6), 
        
        # oceanfloor solute boundary conditions
        ("initial_value", "oceanfloor.O2_conc", 
            [
                0.250,  0.125, 0.050, 0.025, 0.0,   # 5 oxic high SO4
                0.250,  0.125, 0.050, 0.025, 0.0,   # 5 oxic high SO4, no bioturbation
                0.0, 0.0, 0.0, 0.0, 0.0
            ]
        ), # mol m-3
        ("initial_value", "oceanfloor.SO4_conc",  
            [
                28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, # 5 oxic high SO4
                28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, 28756.0e-3, # 5 oxic high SO4, no bioturbation
                10000e-3, 5000e-3, 1000e-3, 100e-3, 0.0
            ]
        ), # mol m-3
             
        # Remineralization
        # sensitivity test for low values for limiting conc
        # ("set_par", "sediment", "reminsed", "FeIIIOxreminlimit", 250.0), # (mol m-3) default = 100e-6*2.5e6 (100 umol Fe(OH)3 g-1 assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3), Van Cappellen 1996)
        # ("set_par", "sediment", "reminsed", "FeIIIOxreminlimit", 0.1*250.0), # (mol m-3)
        ("set_par", "sediment", "reminsed", "FeIIIOxreminlimit", 0.5*250.0), # (mol m-3)

        # Fe-S system
        # Pyrite formation rate
        # ("set_par", "sediment", "PyrH2S", "R_Pyr_H2S", 0.0), # zero rate (disable) pyrite formation
        ("set_par", "sediment", "PyrH2S", "R_Pyr_H2S", 1e2), # 1e5 M yr-1 = 1e5*1e-3 (mol m-3) yr-1, Dale (2015)
        # ("set_par", "sediment", "PyrH2S", "R_Pyr_H2S", 1e1), # sensitivity test for / x10 value
        # Pyrite oxidation rate
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "R_FeS2pyr_O2",  1.0),  # (mol m-3)-1 yr-1,  1e3 M-1 yr-1 = 1e3*1e-3, Dale (2015)
        ("set_par", "sediment", "redox_FeS2pyr_O2", "R_FeS2pyr_O2",  10.0),  # (mol m-3)-1 yr-1  sensitivity test: x10 pyrite oxidation rate
        # ("set_par", "sediment", "redox_FeS2pyr_O2", "R_FeS2pyr_O2",  0.0), # disable pyrite oxidation

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
gr(size=(1500, 1000))
pager = PALEOmodel.PlotPager((3,5), (legend_background_color=nothing, margin=(5, :mm)))

colrange=1:num_columns # number of columns

# check numerics: DAE constraint variables ~0
#  plot_tracers(run.output; colT=[1.0, last(tspan)], tracers=["SmIIaqtot_constraint", "FeIIaqtot_constraint"], colrange, pager=pager)

plot_Corg_O2(run.output; Corgs=["Corg1", "Corg2"], colT=[first(tspan), last(tspan)], colrange, pager)
plot_solutes(run.output; colT=[first(tspan), last(tspan)], solutes=["P", "DIC", "TAlk", "SO4", "SmIIaqtot", "CH4", "H2", "FeIIaqtot"], colrange, pager)
plot_sediment_FeS_summary(run.output; colrange, pager)
plot_solids(run.output; colT=[first(tspan), last(tspan)], solids=["FeHR", "FeMR", "FePR", "FeSm", "FeS2pyr"], colrange, pager)
plot_rates(run.output; colT=[first(tspan), last(tspan)], remin_rates=["reminOrgOxO2", "reminOrgOxFeIIIOx", "reminOrgOxSO4", "reminOrgOxCH4"], colrange, pager)
plot_carbchem(run.output; include_constraint_error=true, colT=last(tspan), colrange, pager)

pager(:newpage)


# summary plots of solute fluxes vs O2_conc
O2_conc = [PALEOmodel.get_array(run.output, "oceanfloor.O2_conc", (tmodel=1e12, cell=i)).values for i in 1:num_columns]
SO4_conc = [PALEOmodel.get_array(run.output, "oceanfloor.SO4_conc", (tmodel=1e12, cell=i)).values for i in 1:num_columns]
soluteflux_O2 = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_O2", (tmodel=1e12, cell=i)).values for i in 1:num_columns]
soluteflux_FeII = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_FeIIaqtot", (tmodel=1e12, cell=i)).values for i in 1:num_columns]
soluteflux_H2S = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_SmIIaqtot", (tmodel=1e12, cell=i)).values for i in 1:num_columns]
soluteflux_CH4 = [PALEOmodel.get_array(run.output, "fluxOceanfloor.soluteflux_CH4", (tmodel=1e12, cell=i)).values for i in 1:num_columns]

reminOrgOxO2_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxO2", (tmodel=1e12, column=i)).values) for i in 1:num_columns]
reminOrgOxFeIIIox_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxFeIIIOx", (tmodel=1e12, column=i)).values) for i in 1:num_columns]
reminOrgOxSO4_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxSO4", (tmodel=1e12, column=i)).values) for i in 1:num_columns]
reminOrgOxCH4_total = [sum(PALEOmodel.get_array(run.output, "sediment.reminOrgOxCH4", (tmodel=1e12, column=i)).values) for i in 1:num_columns]

redox_H2S_O2_total = [sum(PALEOmodel.get_array(run.output, "sediment.redox_H2S_O2", (tmodel=1e12, column=i)).values) for i in 1:num_columns]
redox_CH4_SO4_total = [sum(PALEOmodel.get_array(run.output, "sediment.redox_CH4_SO4", (tmodel=1e12, column=i)).values) for i in 1:num_columns]

# burial fluxes
burial_FeS2pyr = [sum(PALEOmodel.get_array(run.output, "fluxOceanBurial.flux_FeS2pyr", (tmodel=1e12, cell=i)).values) for i in 1:num_columns]



gr(size=(1200, 600))
pager = PALEOmodel.PlotPager((2, 3), (legend_background_color=nothing, margin=(5, :mm)))

# define columns to plot rate summaries
column_ids = [
    #       x axis var      title suffix            x axis label
    (1:5,   O2_conc,        "[SO4] = 28mM",         "[O2] (mol m-3)"),
    (6:10,  O2_conc,        "[SO4] = 28mM, no bio", "[O2] (mol m-3)"),
    (11:15,  SO4_conc,       "[O2] = 0",             "[SO4] (mol m-3)"),
]

for (crange, xvar, titlesuffix, xlabel) in column_ids  
    p = plot(title="solute fluxes, $titlesuffix", xlabel=xlabel, ylabel="flux (mol m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true)
    plot!(p, xvar[crange], soluteflux_O2[crange]; label="O2")
    plot!(p, xvar[crange], soluteflux_FeII[crange]; label="FeII")
    plot!(p, xvar[crange], soluteflux_H2S[crange]; label="H2S")
    plot!(p, xvar[crange], soluteflux_CH4[crange]; label="CH4")
    
    pager(p)
end

for (crange, xvar, titlesuffix, xlabel) in column_ids  
    p = plot(title="remin fluxes, $titlesuffix", xlabel=xlabel, ylabel="flux (mol O2eq m-2 yr-1)", ylim=(-2.0, 1.0), xflip=true);
    plot!(p, xvar[crange], reminOrgOxO2_total[crange]; label="O2")
    plot!(p, xvar[crange], reminOrgOxFeIIIox_total[crange]; label="FeIII")
    plot!(p, xvar[crange], reminOrgOxSO4_total[crange]; label="SO4")
    plot!(p, xvar[crange], reminOrgOxCH4_total[crange]; label="CH4")
    plot!(p, xvar[crange], redox_H2S_O2_total[crange]; label="redox H2S + O2")
    pager(p)
end
   
for (crange, xvar, titlesuffix, xlabel) in column_ids   
    p = plot(title="H2S fluxes, $titlesuffix", xlabel=xlabel, ylabel="flux (mol S m-2 yr-1)", ylim=(-1.0, 1.0), xflip=true);
    plot!(p, xvar[crange], -0.5.*reminOrgOxSO4_total[crange]; label="remin SO4");
    plot!(p, xvar[crange], redox_CH4_SO4_total[crange]; label="CH4 + SO4 -> H2S");
    plot!(p, xvar[crange], -1.0.*soluteflux_H2S[crange]; label="solute H2S");
    plot!(p, xvar[crange], -2.0.*burial_FeS2pyr[crange]; label="pyrite burial");
    pager(p)
end

for (crange, xvar, titlesuffix, xlabel) in column_ids 
    p = plot(title="x = S-II buried / SR, $titlesuffix", xlabel=xlabel, ylabel="x = S-II buried / SR", ylim=(0, 1.0), xflip=true);
    plot!(p,xvar[crange], (2.0.*burial_FeS2pyr[crange]) ./ (-0.5.*reminOrgOxSO4_total[crange]) ; label="x");
    pager(p)
end

pager(:newpage)

