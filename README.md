# PALEOsediment.jl

[![CI](https://github.com/PALEOtoolkit/PALEOsediment.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOsediment.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOsediment.jl/dev)

PALEOtoolkit sediment components and standalone examples

**NB: work-in-progress - this repo contains initial minimal examples only to test infrastructure.**

## Installation and running a minimal example

### Installation

NB: requires Julia 1.9 or later.  To check the Julia version:

    julia> versioninfo()

Clone this github repository to local directory `PALEOsediment`: from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOsediment.jl.git PALEOsediment

Start julia and navigate to the `PALEOsediment/examples` folder, and run `setup.jl` to configure the `PALEOsediment/examples`
Julia environment to use the local (downloaded) version of the PALEOsediment package:

    julia> cd("PALEOsediment/examples")
    julia> include("setup.jl") # use the local version of PALEOsediment packages to allow local modifications
   
### Running a minimal example
Start julia and navigate to the `PALEOsediment` folder, then:

    julia> cd("examples/boudreau1996")
    julia> import Pkg
    julia> Pkg.activate("..") # use the PALEOsediment/examples environment

    julia> include("PALEO_examples_sediment_x10.jl")


## Using PALEOsediment Reactions from other models

The PALEO Reactions comprising the PALEOsediment models are available when the registered PALEOsediment package is loaded (without downloading the repository), ie

    julia> Pkg.add("PALEOsediment")
    julia> import PALEOsediment


