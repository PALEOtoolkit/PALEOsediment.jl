
function plot_Corg_O2(
    output;
    Corgs=["Corg1", "Corg2"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for icol in colrange
        pager(plot(title="Temperature (K) $icol", output, ["sediment.temp"], (tmodel=colT, column=icol), xlabel="Temperature (K)", swap_xy=true))
    end
    for icol in colrange
        pager(plot(title="[O2] $icol", output, ["sediment.O2_conc"], (tmodel=colT, column=icol), xlabel="O2 conc (mol m-3 solute)", swap_xy=true))
    end
    for icol in colrange
        pager(plot(title="[Corg] $icol", output, "sediment.".*Corgs.*"_conc", (tmodel=colT, column=icol), xlabel="Corg conc (mol m-3 solid)", swap_xy=true))
    end
    for icol in colrange
        pager(plot(title="O2 flux Oceanfloor $icol", output, ["fluxOceanfloor.soluteflux_O2"], (cell=icol,), ylabel="O2 (mol yr-1)",))
    end
    for Corg in Corgs
        for icol in colrange
            pager(plot(title="$Corg input output $icol", output, ["fluxOceanfloor.particulateflux_$Corg", "fluxOceanBurial.flux_$Corg"], 
                (cell=icol,), ylabel="Corg (mol yr-1)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_Corg_RCmultiG(
    output;
    Corg_indices=1:14,
    colT=1e12,
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    Corg_vars = vcat(
        ["Corg_conc"],
        ["Corg_$(i)_conc" for i in Corg_indices],
    )
   
    for icol in colrange
        pager(plot(title="[Corg] $icol", output, "sediment.".*Corg_vars, (tmodel=colT, column=icol), labellist=copy(Corg_vars), xlabel="Corg conc (mol m-3 solid)", swap_xy=true))
    end

    pager(:newpage)

    return nothing
end

function plot_phi(
    output;
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        pager(plot(title="porosity $icol", output, ["sediment.phi_solid", "sediment.phi"], (column=icol, tmodel=1e12), xlabel="volume fraction", swap_xy=true))
    end

end

function plot_w(
    output;
    colT=[1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        pager(plot(title="advection velocity $icol", output, ["sediment.w_solid_upper", "sediment.w_solid_lower", "sediment.w_solute_upper",  "sediment.w_solute_lower"], (tmodel=colT, column=icol,), xlabel="w (m yr-1)", swap_xy=true))
    end

end

function plot_biorates(
    output;
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        pager(plot(title="bioturbation $icol", output, "sediment.diff_bioturb", (column=icol, tmodel=1e12), label="sediment.diff_bioturb", xlabel="effective diffusivity (m^2 yr-1)", swap_xy=true))
    end

    for icol in colrange
        pager(plot(title="bioirrigation $icol", output, "sediment.alpha_bioirrig", (column=icol, tmodel=1e12), label="sediment.alpha_bioirrig", xlabel="solute exchange (yr-1)", swap_xy=true))
    end

end


function plot_solutes(
    output;
    solutes=["P", "SO4", "H2S", "CH4"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solutes
        for icol in colrange
            pager(plot(title="[$s] $icol", output, ["sediment.$(s)_conc"], (tmodel=colT, column=icol), xlabel="$s conc (mol m-3 solute)", swap_xy=true))
        end
    end

    for s in solutes
        for icol in colrange
            pager(plot(title="$s flux Oceanfloor $icol", output, ["fluxOceanfloor.soluteflux_$s"], (cell=icol,), ylabel="$s (mol yr-1)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_conc_summary(
    output;
    species=["P", "SO4", "H2S", "CH4"],
    xlims=(0, Inf),
    xscale=:identity,
    colT=1e12,
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    
    species_vars = "sediment.".*species.*"_conc"
    for icol in colrange
        pager(
            plot(
                title="concentrations $icol t=$colT (yr)", output, species_vars, (tmodel=colT, column=icol); 
                legend=:bottomleft, labellist=copy(species), xlabel="conc (mol m-3)", swap_xy=true, xscale, xlims,
            )
        )
    end

    pager(:newpage)

    return nothing
end

function plot_solute_deltas(
    output;
    solutes=["SO4", "H2S"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solutes
        for icol in colrange
            pager(plot(title="[$s] delta $icol", output, ["sediment.$(s)_conc.v_delta"], (tmodel=colT, column=icol), xlabel="$s delta (per mil)", swap_xy=true))
        end
    end

    for s in solutes
        for icol in colrange
            pager(plot(title="$s flux Oceanfloor delta $icol", output, ["fluxOceanfloor.soluteflux_$s.v_delta"], (cell=icol,), ylabel="$s delta (per mil)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_solid_deltas(
    output;
    solids=["Corg"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solids
        for icol in colrange
            pager(plot(title="[$s] delta $icol", output, ["sediment.$(s)_conc.v_delta"], (tmodel=colT, column=icol), xlabel="$s delta (per mil)", swap_xy=true))
        end
    end

    for s in solids
        for icol in colrange
            pager(plot(title="$s flux Oceanfloor delta $icol", output, ["fluxOceanfloor.particulateflux_$s.v_delta"], (cell=icol,), ylabel="$s delta (per mil)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_solids(
    output;
    solids=["Corg1"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solids
        for icol in colrange
            pager(plot(title="[$s] $icol", output, ["sediment.$(s)_conc"], (tmodel=colT, column=icol), xlabel="$s conc (mol m-3 solid phase)", swap_xy=true))
        end
    end

    for s in solids
        for icol in colrange
            pager(plot(title="$s flux Oceanfloor $icol", output, ["fluxOceanfloor.particulateflux_$s"], (cell=icol,), ylabel="$s (mol yr-1)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_solids_volume_frac(
    output;
    solids=["Corg1"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solids
        for icol in colrange
            pager(plot(title="[$s] $icol", output, ["sediment.$(s)_volume_frac"], (tmodel=colT, column=icol), xlabel="solid phase volume fraction", swap_xy=true))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_tracers(
    output;
    tracers=["SmIIaqtot_constraint", "FeIIaqtot_constraint"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for st in tracers
        for icol in colrange
            pager(plot(title="[$st] $icol", output, ["sediment.$(st)"], (tmodel=colT, column=icol), xlabel="$st", swap_xy=true))
        end
    end

    pager(:newpage)

    return nothing
end

"""
    normalize_by_volume(rate_var, volume_var) -> rate_var / volume_var

Convert per-cell rate (mol yr-1) to `rate_var / volume_var` (mol m-3 yr-1)
"""
function normalize_by_volume(rate_var::PALEOmodel.FieldArray, volume_var::PALEOmodel.FieldArray)
    return PALEOmodel.FieldArray(rate_var.name*" / volume", rate_var.values./volume_var.values, rate_var.dims, rate_var.attributes)
end


function plot_rates(
    output;
    remin_rates=["reminOrgOxO2", "reminOrgOxSO4", "reminOrgOxCH4"],
    plot_freminOrgTot=false,
    colT=1e12,
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        p = plot()
        # remin_Corg at each time
        for t in colT
            volume = PALEOmodel.get_array(output, "sediment.volume", (tmodel=t, column=icol))
            remin_Corg = PALEOmodel.get_array(output, "sediment.remin_Corg", (tmodel=t, column=icol))
            remin_Corg_norm = normalize_by_volume(remin_Corg, volume)
            plot!(p, remin_Corg_norm; swap_xy=true)
        end

        plot!(p, title="Corg remin rate $icol", xlabel="remin (mol Corg m-3 yr-1)")

        pager(p)
    end

    if plot_freminOrgTot
        for icol in colrange
            pager(plot(output, "sediment.freminOrgTot", (tmodel=colT, column=icol); title="ox dep rate factor $icol", xlim=(0, Inf), swap_xy=true))
        end
    end

    for icol in colrange
        # remin_rates at last time
        volume = PALEOmodel.get_array(output, "sediment.volume", (tmodel=last(colT), column=icol))
        p = plot()

        for r in remin_rates
            rate = PALEOmodel.get_array(output, "sediment.".*r, (tmodel=last(colT), column=icol))
            rate_norm = normalize_by_volume(rate, volume)
            plot!(p, rate_norm; swap_xy=true)
        end
        plot!(p; title="Org matter ox rate $icol", xlabel="ox rate (mol O2eq m-3 yr-1)")

        pager(p)
    end

end

function plot_sediment_FeS_summary(
    output;
    SmII_species = ["TH2S", "H2S", "HS", "FeSaq",],
    FeII_species = ["FeIIsolutetot", "FeII", "FeSaq", "FeIIadsorb", "FeIItot"],
    colT=last(PALEOmodel.get_array(output, "global.tmodel").values), # model time for column plots
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)

    for icol in colrange
        pager(
            plot(title="S-II speciation $icol t=$colT (yr)", output, "sediment.".*SmII_species.*"_conc", 
                (tmodel=colT, column=icol);
                swap_xy=true, xlabel="species conc (mol m-3)", labellist=copy(SmII_species), plotargs...)
            )
    end

    for icol in colrange
        pager(
            plot(title="FeII speciation $icol t=$colT (yr)", output, "sediment.".*FeII_species.*"_conc", 
                (tmodel=colT, column=icol);
                swap_xy=true, xlabel="species conc (mol m-3)", labellist=copy(FeII_species), xscale=:log10, plotargs...)
            )
    end

    for icol in colrange
        pager(
            plot(title="FeS saturation $icol t=$colT (yr)", output, "sediment.Omega_FeSaq", 
                (tmodel=colT, column=icol);
                swap_xy=true, xlabel="Omega_FeSaq", labellist=["Omega_FeSaq"], plotargs...)
            )
    end
   
    return nothing
end

function plot_carbchem(
    output;
    include_constraint_error=false,
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    
    for icol in colrange
        pager(plot(title="pHtot $icol", output, ["sediment.pHtot"], (tmodel=colT, column=icol), xlabel="pH (total)", swap_xy=true))
    end

    if include_constraint_error
        for icol in colrange
            # DAE constraint (should be ~zero, difference between required and calculated TAlk)
            pager(plot(title="TAlk_constraint $icol", output, ["sediment.TAlk_constraint"], (tmodel=colT, column=icol), xlabel="mol", swap_xy=true))
        end
    end

    for icol in colrange
        pager(plot(title="Omega Ca carbonate $icol", output, ["sediment.OmegaCA", "sediment.OmegaAR"], (tmodel=colT, column=icol), xlabel="Omega", swap_xy=true))
    end
    
    pager(:newpage)

    return nothing
end

function plot_FeP(
    output;
    FeP=["PFeHR_theta", "PFeMR_theta"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
   
    for fp in FeP
        for icol in colrange
            pager(plot(title="[$fp] $icol", output, ["sediment.$(fp)"], (tmodel=colT, column=icol), xlabel="$fp P:Fe (mol/mol)", swap_xy=true))
        end
    end

    pager(:newpage)

    return nothing
end


function plot_budget(
    output;
    name="Mn",
    solids=["MnHR", "MnMR"],
    solutes=["MnII"],
    extras=[],
    solute_burial_flux=false,
    stoich_factors=Dict(), # eg Dict("Corg"=>1/106) for P:Corg ratio
    concxscale=:log10, 
    concxlims=(1e-3, Inf),
    fluxyscale=:identity, 
    fluxylims=(-Inf, Inf),
    colT=last(PALEOmodel.get_array(output, "global.tmodel").values), # model time for column plots
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    # scale by stoich_factors[species] 
    function scale_stoich(Cconc, species, label)
        if haskey(stoich_factors, species)
            sf = stoich_factors[species]
            Cconc *= sf
            label *= " * $(round(sf, sigdigits=3))"
        end
        return Cconc, label
    end

    # plot concentrations
    all_species=vcat(solutes, solids)
    for icol in colrange
        p = plot(title="$name concentrations $icol")        
        for sp in all_species
            conc = PALEOmodel.get_array(output, "sediment.$(sp)_conc", (column=icol, tmodel=colT))
            conc, label = scale_stoich(conc, sp, sp)
            plot!(p, conc; label, swap_xy=true)
        end
        plot!(p, xlabel="conc (mol m-3)", xscale=concxscale, xlims=concxlims)
        pager(p)
    end

    # plot input/output fluxes and total budget 
    for icol in colrange
        tmodel = PALEOmodel.get_array(output, "sediment.tmodel").values
        total_flux_into_sediment = zeros(length(tmodel))
        p = plot(title="$name budget $icol")
        for s in solutes
            # +ve is sediment -> ocean
            fluxname = "fluxOceanfloor.soluteflux_$s"
            solute_flux = PALEOmodel.get_array(output, fluxname, (cell=icol,))
            solute_flux, label = scale_stoich(solute_flux, s, "-"*fluxname)
            total_flux_into_sediment .-= solute_flux.values
            plot!(p, -1*solute_flux; label, linestyle=:dash)
        end
        if solute_burial_flux
            for s in solutes
                # +ve is out of sediment
                fluxname = "fluxOceanBurial.flux_$s"
                burial_flux = PALEOmodel.get_array(output, fluxname, (cell=icol,))
                burial_flux, label = scale_stoich(burial_flux, s, "-"*fluxname)
                total_flux_into_sediment .-= burial_flux.values
                plot!(p, -1*burial_flux; label)
            end
        end
        for s in solids
            # +ve is into sediment
            fluxname = "fluxOceanfloor.particulateflux_$s"
            particulate_flux = PALEOmodel.get_array(output, fluxname, (cell=icol,))
            particulate_flux, label = scale_stoich(particulate_flux, s, fluxname)
            total_flux_into_sediment .+= particulate_flux.values
            plot!(p, particulate_flux; label)
        end
        for s in solids
            # +ve is out of sediment
            fluxname = "fluxOceanBurial.flux_$s"
            burial_flux = PALEOmodel.get_array(output, fluxname, (cell=icol,))
            burial_flux, label = scale_stoich(burial_flux, s, "-"*fluxname)
            total_flux_into_sediment .-= burial_flux.values
            plot!(p, -1*burial_flux; label)
        end

        for x in extras
            # +ve is in to sediment
            x_flux = PALEOmodel.get_array(output, x, (cell=icol,))
            x_flux, label = scale_stoich(x_flux, x, x)
            total_flux_into_sediment .+= x_flux.values
            plot!(p, x_flux; label)
        end

        plot!(p, tmodel, total_flux_into_sediment, label="total", color=:black)

        plot!(p, ylabel="$name flux (mol yr-1)", yscale=fluxyscale, ylims=fluxylims)
        pager(p)
    end
    
    pager(:newpage)

    return nothing
end

function plot_summary_stacked(
    output;
    species=["P", "SO4", "H2S", "CH4"],
    others=[],
    budgets=[],
    plot_pHfree=true,
    xscale=:identity,
    xlims=(0, Inf),
    pHxlims=(-Inf, Inf),
    oxscale=:identity,
    oxlims=(0, Inf),
    budgetylims=(-1e-6, 1e-6),    
    ylims=(-Inf, Inf),
    colT=colT=last(PALEOmodel.get_array(output, "global.tmodel").values), # model time for column plots
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in species
        pager(
            plot(
                title="$s t=$colT (yr)", output, "sediment.".*s.*"_conc", (tmodel=colT, column=collect(colrange)); 
                legend=:bottomleft, labellist=string.(collect(colrange)), xlabel="$s conc (mol m-3)", swap_xy=true, xscale, xlims, ylims,
            )
        )
    end

    for o in others
        pager(
            plot(
                title="$o t=$colT (yr)", output, "sediment.".*o, (tmodel=colT, column=collect(colrange)); 
                legend=:bottomleft, labellist=string.(collect(colrange)), swap_xy=true, oxscale, oxlims, ylims,
            )
        )
    end

    if plot_pHfree
        pager(
            plot(
                title="pHfree t=$colT (yr)", output, "sediment.pHfree", (tmodel=colT, column=collect(colrange)); 
                legend=:bottomleft, labellist=string.(collect(colrange)), xlabel="pH (free)", swap_xy=true, xlims=pHxlims, ylims,
            )
        )
    end

    if !isempty(budgets)
        for i in colrange
            pager(
                plot(title="budget check $i", output, budgets, (cell=i,), ylabel="net input (mol yr-1)", ylims=budgetylims)
            )
        end
    end

    return nothing
end

function plot_budgets(
    output;
    budgets=[],
    ylims=(-1e-6, 1e-6),
    clamp_y_fac = 100.0, # workaround for Plots.jl bug: values too much larger than ylims result in whole line being omitted
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
   
    for i in colrange
        pager(
            plot(title="budget check $i", output, budgets, (cell=i,); 
                ylabel="net input (mol yr-1)", 
                ylims, map_values=y->clamp(y, clamp_y_fac*ylims[1], clamp_y_fac*ylims[2])),
        )
    end

    return nothing
end