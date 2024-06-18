# PALEOsediment.jl

[![CI](https://github.com/PALEOtoolkit/PALEOsediment.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOsediment.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOsediment.jl/dev)

Sediment components and standalone examples for the [PALEOtoolkit](https://github.com/PALEOtoolkit) biogeochemical framework.

The `PALEOsediment` package provides an implementation of sediment transport (PALEO reaction `ReactionSedimentTransport`) for n x 1D sediment columns.

This can then be combined with biogeochemistry implemented by the `PALEOaqchem` and `PALEOboxes` packages and solvers from the `PALEOmodel` package  to create either standalone sediment models or coupled water-column - sediment configurations.

The implementation of 1D sediment reaction-transport is standard (eg [Van Cappellen & Wang (1996)](https://dx.doi.org/10.2475/ajs.296.3.197), [Boudreau (1996)](https://dx.doi.org/10.1016/0098-3004(95)00115-8)) and uses the PALEOtoolkit framework to provide:
- extensible configuration for biogeochemical species and reactions defined in a .yaml configuration file
- standalone and coupled sediment configurations
- fully implicit numerical solution using a combination of sparse automatic differentation to generate Jacobians, and efficient solvers including pseudo-transient-continuation for steady-state solutions.

## Documentation

Documentation is available online at https://paleotoolkit.github.io/PALEOsediment.jl/

## Installation

### Using PALEOsediment Reactions from other models

The PALEOsediment Reactions are available to [PALEOtoolkit](https://github.com/PALEOtoolkit) when the registered PALEOsediment package is installed and loaded:

    julia> Pkg.add("PALEOsediment") # install PALEOsediment package in the currently active Julia environment
    julia> import PALEOsediment

### Running the PALEOsediment examples

To install and run the PALEOsediment examples, clone this github repository to local directory `PALEOsediment` and run the examples from the Julia REPL.

Quickstart assuming a recent Julia installation: from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOsediment.jl.git PALEOsediment

Start julia and navigate to the `PALEOsediment/examples` folder, and run `setup.jl` to configure the `PALEOsediment/examples`
Julia environment to use the local (downloaded) version of the PALEOsediment package:

    julia> cd("PALEOsediment/examples")
    julia> include("setup.jl") # use the local version of PALEOsediment packages to allow local modifications
   
Examples are in subfolders of `PALEOsediment/examples` and the use the `PALEOsediment/examples` Julia environment.


See [Installation and getting started](https://paleotoolkit.github.io/PALEOtutorials.jl/dev/ExampleInstallConfig/)
in the [PALEOtutorials](https://github.com/PALEOtoolkit/PALEOtutorials.jl) repository for more details including installation and configuration of Julia.
