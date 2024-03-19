using Test
using Logging
import DataFrames
import NLsolve
import LineSearches

import PALEOboxes as PB

import PALEOsediment
import PALEOaqchem
import PALEOmodel


@testset "sediment boudreau1996" begin

skipped_testsets = [
    # "Corg_O2",    
]

!("Corg_O2" in skipped_testsets) && @testset "Corg_O2" begin

    include("config_sediment_expts.jl")

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_sediment_cfg.yaml"), 
        "sediment_Corg_O2",
    )
        
    config_sediment_expts(model, 
        [
            "baseline",
        ]
    ) 

    tspan=(0,10000.0)
    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # Newton, no(?) line search
    sol = PALEOmodel.SteadyState.steadystateForwardDiff(
        run, initial_state, modeldata, 0.0;
        solvekwargs=(
            ftol=10e-9,
            # ftol=5e-5,
            method=:newton, 
            store_trace=true, 
            show_trace=true
        ),
    );
    
    @test NLsolve.converged(sol)

    println("conservation checks:")
    conschecks = [ ]    
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(run.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    
    # check values with quadratic grid, 1000 bins
    checkvals = [   
        ("fluxOceanfloor", "soluteflux_O2",  [-1.64229, -0.24071, -0.23907],      1e-3),
        ("fluxOceanfloor", "soluteflux_H2S", [0.025073, 6.8626e-5, 0.00088524],   2.5e-2),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end



end
