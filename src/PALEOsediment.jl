module PALEOsediment

import PALEOboxes as PB

function moduledir()
    return dirname(@__DIR__)
end

include("sediment/Sediment.jl")

end # module PALEOsediment
