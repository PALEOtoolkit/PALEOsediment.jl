
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

function plot_w(
    output;
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)

    for icol in colrange
        pager(plot(title="advection velocity $icol", output, ["sediment.w_solid", "sediment.w_solute"], (column=icol, tmodel=1e12), xlabel="w (m yr-1)", swap_xy=true))
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

function plot_solute_deltas(
    output;
    solutes=["SO4", "H2S"],
    colT=[-1e12, 0.1, 1.0, 10.0, 100.0, 1000.0, 1e12],
    colrange=1:3,
    pager=PALEOmodel.DefaultPlotPager(),
)
    for s in solutes
        for icol in colrange
            pager(plot(title="[$s] delta $icol", output, ["sediment.$(s).v_delta"], (tmodel=colT, column=icol), xlabel="$s delta (per mil)", swap_xy=true))
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


function plot_rates(
    output;
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
                ["sediment.reminOrgOxO2", "sediment.reminOrgOxSO4", "sediment.reminOrgOxCH4"], 
                (tmodel=colT[end], column=icol), # only last time
                xlabel="ox rate (mol O2eq yr-1)",
                swap_xy=true,
            )
        )
    end

end
