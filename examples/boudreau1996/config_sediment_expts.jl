
"test cases and examples for sediment"
function config_sediment_expts(
    model, expts; 
)

    for expt in expts        
        println("Add expt: ", expt)
        if expt == "baseline"
            # defaults

        elseif length(expt) == 5 && expt[1] == "set_par"
            # generic parameter set (set_par, <domain>, <reaction>, <parname>, <parvalue)
            _, domname, reactname, parname, parvalue = expt            
            PB.set_parameter_value!(model, domname, reactname, parname, parvalue)

        elseif length(expt) == 3 && expt[1] == "initial_value"
            # generic set variable initial value
            _, domvarname, varval = expt           
            PB.set_variable_attribute!(model, domvarname, :initial_value, varval)
            
        elseif expt == "slowredox"
            # reduce rate of secondary redox reactions to help convergence
            dom_sediment = PB.get_domain(model, "sediment")
            redox_H2S_O2 = PB.get_reaction(dom_sediment, "redox_H2S_O2")
            PB.setvalue!(redox_H2S_O2.pars.R_H2S_O2, 1.8e3)
            @info "expt slowredox set redox_H2S_O2.pars.R_H2S_O2 = $(redox_H2S_O2.pars.R_H2S_O2[])"
            redox_CH4_O2 = PB.get_reaction(dom_sediment, "redox_CH4_O2")
            PB.setvalue!(redox_CH4_O2.pars.R_CH4_O2, 1e4)
            @info "expt slowredox set redox_CH4_O2.pars.R_CH4_O2 = $(redox_CH4_O2.pars.R_CH4_O2[])"
        else
            error("unrecognized expt='$(expt)'")
        end
    end

    return nothing
end

