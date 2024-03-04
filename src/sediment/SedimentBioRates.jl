module SedimentBioRates

import SpecialFunctions

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionSedimentBioRates

Calculate sediment bioturbation and bioirrigation rates (for use by sediment transport).

Literature compilation of functional forms for dependency on oceanfloor oxygen, Corg flux, and depth within sediment, from
[Boudreau1996](@cite), [Archer2002](@cite), [Arndt2011](@cite), [Dale2015](@cite)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSedimentBioRates{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("f_bioTurbRate", "Prescribed", 
            allowed_values=["Prescribed", "CorgArcher2002"],
            description="functional form for bioturbation max rate"),

        PB.ParString("f_bioTurbDepth", "ConstCutoff",
            allowed_values=["ConstCutoff", "ErfcCutoff", "Exp1Cutoff", "Exp2Cutoff"],
            description="functional form of bioturbation rate with depth in sediment"),

        PB.ParString("f_bioIrrigDepth", "ConstCutoff",
            allowed_values=["ConstCutoff", "ErfcCutoff", "Exp1Cutoff", "Exp2Cutoff"],
            description="functional form of bioirrigation rate with depth in sediment"),

        PB.ParString("f_bioO2", "None",
            allowed_values=["None", "MM", "Dale2015"],
            description="functional form of bioturbation and bioirrigation sensitivity to oceanfloor oxygen"),
        PB.ParDouble("bioO2halfmax", 20e-3, units="mol m-3",
            description="oceanfloor [O2] for 50% decrease in bioturbation/bioirrigation"),
        PB.ParDouble("bioO2decreaserate", 12e-3, units="mol m-3",
            description="oceanfloor [O2] sharpness of decrease in bioturbation/bioirrigation"),
        
        PB.ParBool("separate_zbio", false,
            description="true to use zbioturb, zbioirrig variables, false to use zbio for both bioturbation and bioirrigation"),
    )


end

function PB.register_methods!(rj::ReactionSedimentBioRates, model::PB.Model)

    vars = [
        PB.VarDepColumn("oceanfloor_Dbio"=>"oceanfloor.Dbio",
            "m^2 yr-1", "characteristic value of bioturbation effective diffusivity"),
        PB.VarDepColumn("oceanfloor_O2_conc"=>"(oceanfloor.O2_conc)",
            "mol m-3", "oceanfloor oxygen concentration"),
        PB.VarDepColumn("oceanfloor_alpha"=>"oceanfloor.alpha",
            "yr-1", "characteristic value for bioirrigation exchange rate"),
        PB.VarDep("zmid",
            "m",        "mean depth of box"),
        PB.VarProp("diff_bioturb",
            "m^2 yr-1", "bioturbation effective diffusivity"),
        PB.VarProp("alpha_bioirrig",
            "yr-1",     "bioirrigation exchange rate"),
    ]
    if rj.pars.separate_zbio[]
        append!(vars, [
            PB.VarDepColumn("oceanfloor_zbioturb"=>"oceanfloor.zbioturb",
                "m", "characteristic depth for bioturbation (+ve)"),
            PB.VarDepColumn("oceanfloor_zbioirrig"=>"oceanfloor.zbioirrig",
                "m", "characteristic depth for bioirrigation (+ve)"),
        ])
    else
        append!(vars, [
            PB.VarDepColumn("oceanfloor_zbio"=>"oceanfloor.zbio",
                "m", "characteristic depth for bioturbation and bioirrigation (+ve)"),
        ])
    end

    funcs = (
        f_bioTurbDepth  =   _get_f_bioDepth(rj.pars.f_bioTurbDepth[]),
        f_bioO2         =   _get_f_bioO2(rj.pars.f_bioO2[]),
        f_bioIrrigDepth =   _get_f_bioDepth(rj.pars.f_bioIrrigDepth[]),
    )

    
    PB.setfrozen!(rj.pars.f_bioTurbDepth)
    PB.setfrozen!(rj.pars.f_bioO2)
    PB.setfrozen!(rj.pars.f_bioIrrigDepth)

    PB.add_method_do!(
        rj,
        do_sediment_bio_rates,
        (PB.VarList_namedtuple(vars),),
        p=funcs,
    )

    if rj.pars.f_bioTurbRate[] == "CorgArcher2002"
        # NB: Variables and method in oceanfloor Domain
        btvars = [
            PB.VarProp("Dbio", 
                "m^2 yr-1", "characteristic value of bioturbation effective diffusivity"),
            PB.VarDep("fluxOceanfloor_Corg"=>"fluxOceanfloor.particulateflux_Corg",
                "mol yr-1", "total oceanfloor organic carbon flux"),
            PB.VarDep("Afloor",
                "m^2", "horizontal area of seafloor at sediment surface"),
        ]
        PB.add_method_do!(
            rj,
            do_sediment_bioturb_Archer,
            (PB.VarList_namedtuple(btvars),);
            domain=PB.get_domain(model, "oceanfloor"),
        )

    end
    PB.setfrozen!(rj.pars.f_bioTurbRate)

    return nothing
end


function _get_f_bioDepth(f_bioDepth)
    f_bioDepthDict = Dict(
        "ConstCutoff" => (zmid, zbio) -> -zmid < zbio,
        # Arndt etal (2011) Eq(56) (a guess, typo in paper?)
        "ErfcCutoff"  => (zmid, zbio) -> 0.5 * SpecialFunctions.erfc(-zmid/zbio - 1.0),
        # Boudreau (1996) Eq(99), the 'second biodiffusion function'
        # Archer etal (2002) Eq(14)
        "Exp2Cutoff"   => (zmid, zbio) -> exp(-zmid^2/(2.0*zbio^2)),
        # Dale (2015) (used for bioirrigation)
        "Exp1Cutoff"   => (zmid, zbio) -> exp(zmid/zbio), # NB: zmid is -ve here
    )
    haskey(f_bioDepthDict, f_bioDepth) || error("unknown f_bioDepth $f_bioDepth")    
    return f_bioDepthDict[f_bioDepth]    
end

function _get_f_bioO2(f_bioO2)
    f_bioO2Dict = Dict(
        "None"       => (pars, vars, icol) -> 1.0,
        # Archer etal (2002) Eq. 14
        "MM"         => (pars, vars, icol) 
                        -> max(vars.oceanfloor_O2_conc[icol], 0.0)/(max(vars.oceanfloorO2_conc[icol], 0.0) + pars.bioO2halfmax[]),
        # Dale etal (2015)
        "Dale2015"   => (pars, vars, icol) 
                        -> 0.5 + 0.5*SpecialFunctions.erf((vars.oceanfloor_O2_conc[icol] - pars.bioO2halfmax[])/pars.bioO2decreaserate[]),
    )
    haskey(f_bioO2Dict, f_bioO2) || error("unknown f_bioO2 $f_bioO2")    
    return f_bioO2Dict[f_bioO2]    
end

function do_sediment_bioturb_Archer(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    
    function f_bioRate_CorgArcher2002(pars, vars, i)
        # Archer (2002) Eq(13)  bioturbation rate as a function of Corg rain rate
        # umol/cm^2/yr                    mol/m^2/yr
        rainCorg = max(vars.fluxOceanfloor_Corg[i]/vars.Afloor[i]*1e6/1e4, 0.0)
        return 0.0232e-4 * pow(rainCorg, 0.85)
    end
    
    for i in cellrange.indices
        vars.Dbio[i] = f_bioRate_CorgArcher2002(pars, vars, i)
    end            

    return nothing
end


function do_sediment_bio_rates(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    
    funcs = m.p

    if pars.separate_zbio[]
        var_zbioturb = vars.oceanfloor_zbioturb
        var_zbioirrig = vars.oceanfloor_zbioirrig
    else
        var_zbioturb = vars.oceanfloor_zbio
        var_zbioirrig = vars.oceanfloor_zbio
    end
    
    for (icol, colindices) in cellrange.columns
        bioTurbRate     = vars.oceanfloor_Dbio[icol] # max rate
        bioIrrigRate    = vars.oceanfloor_alpha[icol]
        zbioturb = var_zbioturb[icol]
        zbioirrig = var_zbioirrig[icol]
        bioO2           = funcs.f_bioO2(pars, vars, icol)

        @inbounds for i in colindices

            bioTurbDepth    = funcs.f_bioTurbDepth(vars.zmid[i], zbioturb)

            vars.diff_bioturb[i] = bioTurbRate * bioTurbDepth * bioO2

            bioIrrigDepth   = funcs.f_bioIrrigDepth(vars.zmid[i], zbioirrig)

            vars.alpha_bioirrig[i] = bioIrrigRate * bioIrrigDepth * bioO2

        end
    end

    return nothing
end

end # module
