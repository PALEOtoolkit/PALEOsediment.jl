module SumColumns_dev

import PALEOboxes as PB

using PALEOboxes.DocStrings


"""
    ReactionSumColumns_dev
 
Sum a quantity `X` by columns, output to `X_sum`

Target variables `X_sum` needs to exist in the boundary domain,
eg created by `ReactionFluxTarget`

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSumColumns_dev{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        
        PB.ParDouble("mult", 1.0, units="",
            description="multiplier to apply X_sum = sum(mult * X)"),
    )

end


function PB.register_methods!(rj::ReactionSumColumns_dev)
    
    @info "ReactionSumColumns_dev.register_methods! $(PB.fullname(rj))"

    vars = [
        PB.VarDep(          "X",     "",    "per-cell quantity to sum over columns"),
        PB.VarContribColumn("X_sum", "",    "column sum of X")
    ]

    PB.add_method_do!(
        rj,
        do_sum_columns,
        (PB.VarList_namedtuple(vars), ),
    )

    return nothing
end

function do_sum_columns(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    mult = pars.mult[]
    for (icol, colindices) in cellrange.columns
        for i in colindices
            vars.X_sum[icol] += mult*vars.X[i]
        end
    end

    return nothing
end



end
