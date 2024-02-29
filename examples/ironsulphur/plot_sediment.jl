
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

function plot_sediment_FeS_summary(
    output;
    colT=last(PALEOmodel.get_array(output, "global.tmodel").values), # model time for column plots
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)

    SmII_species = ["SmIIaqtot", "H2Ssp", "HSm", "FeSaq",]
    for icol in colrange
        pager(
            plot(title="S-II speciation $icol t=$colT (yr)", output, "sediment.".*SmII_species.*"_conc", 
                (tmodel=colT, column=icol);
                swap_xy=true, xlabel="species conc (mol m-3)", labellist=copy(SmII_species), plotargs...)
            )
    end

    FeII_species = ["FeIIaqtot", "FeII", "FeSaq",]
    for icol in colrange
        pager(
            plot(title="FeII speciation $icol t=$colT (yr)", output, "sediment.".*FeII_species.*"_conc", 
                (tmodel=colT, column=icol);
                swap_xy=true, xlabel="species conc (mol m-3)", labellist=copy(FeII_species), plotargs...)
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


function plot_solids(
    output;
    solids=["Fe2O3x"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
   
    for s in solids
        for icol in colrange
            pager(plot(title="[$s] $icol", output, ["sediment.$(s)_conc"], (tmodel=colT, column=icol), xlabel="$s conc (mol m-3 solid)", swap_xy=true))
        end
    end
   
    for s in solids
        for icol in colrange
            pager(plot(title="$s input output $icol", output, ["fluxOceanfloor.particulateflux_$s", "fluxOceanBurial.flux_$s"], 
                (cell=icol,), ylabel="$s (mol yr-1)",))
        end
    end

    pager(:newpage)

    return nothing
end

function plot_rates(
    output;
    remin_rates=["reminOrgOxO2", "reminOrgOxSO4", "reminOrgOxCH4"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        pager(
            plot(
                title="Corg remin rate $icol",
                output,
                ["sediment.remin_Corg"],
                (tmodel=colT, column=icol),
                xlabel="remin (mol Corg yr-1)",
                swap_xy=true
            )
        )
    end

    for icol in colrange
        pager(
            plot(
                title="Org matter ox rate $icol",
                output, 
                "sediment.".*remin_rates, 
                (tmodel=colT[end], column=icol), # only last time
                xlabel="ox rate (mol O2eq yr-1)",
                swap_xy=true,
            )
        )
    end

end
