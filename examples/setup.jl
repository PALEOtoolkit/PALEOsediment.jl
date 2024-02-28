# Setup PALEOsediment/examples environment to use local downloaded PALEOsediment package code
# This only needs to be run once, after cloning the github repository

import Pkg

Pkg.activate(".") # use the PALEOsediment/examples environment
Pkg.develop(path="../")   # use the local version of PALEOsediment packages to allow local modifications
Pkg.add(url="https://github.com/PALEOtoolkit/NLsolve.jl", rev="project_region")
Pkg.instantiate()
Pkg.update()
