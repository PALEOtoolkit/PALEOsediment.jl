using Test
using Documenter


import Pkg

Pkg.add(url="https://github.com/PALEOtoolkit/NLsolve.jl", rev="project_region")

ENV["GKSwstype"] = "100" # to run Plots.jl GR on a headless system

@testset "PALEOsediment all" begin

@testset "sediment examples" begin

    include("../examples/boudreau1996/runtests.jl")

end

doctest(PALEOsediment; manual=false)  

end

delete!(ENV, "GKSwstype"); # undo workaround for Plots.jl on a headless system