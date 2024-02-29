var documenterSearchIndex = {"docs":
[{"location":"collated_examples/boudreau1996/README/#Sediment-Boudreau-(1996)-test-cases","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"","category":"section"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"To run sediment test cases from (Boudreau, 1996), assuming the Julia environment PALEOsediment/examples has already been activated:","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"julia> cd(\"examples/boundreau1996\")\njulia> include(\"PALEO_examples_sediment.jl\")","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"This will run and plot output (NB: the first run will be slow as Julia JIT compiles the code).","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"The configuration includes three sediment columns, with:","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"Shelf / slope case\nRise case\nRise, no bioturbation","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"NB: this configuration doesn't include N, so results differ in detail from the paper.","category":"page"},{"location":"collated_examples/boudreau1996/README/#Oceanfloor-[O2]-and-[SO4]-gradients","page":"Sediment Boudreau (1996) test cases","title":"Oceanfloor [O2] and [SO4] gradients","text":"","category":"section"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"julia> include(\"PALEO_examples_sediment_x10.jl\")","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"The configuration includes ten sediment columns, with the physical environment for the (Boudreau, 1996) shelf/slope case, no bioturbation, and Corg input with two reactivity fractions (a combination of the high-reactivity fraction from the shelf/slope case, and the low reactivity fraction from the rise case).","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"Columns 1-5: oceanfloor [O2] gradient, constant [SO4] = 28mM\nColumns 6-10: oceanfloor [SO4] gradient, at constant [O2] = 0","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"Summary plots show oceanfloor solute fluxes and remineralization pathways:","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"(Image: O2 and SO4 gradient summary figure)","category":"page"},{"location":"collated_examples/boudreau1996/README/#Figure-1","page":"Sediment Boudreau (1996) test cases","title":"Figure 1","text":"","category":"section"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"Oceanfloor solute fluxes and remineralization pathways vs oceanfloor [O2] and [SO4] concentration","category":"page"},{"location":"collated_examples/boudreau1996/README/#Sulphur-isotope-example","page":"Sediment Boudreau (1996) test cases","title":"Sulphur isotope example","text":"","category":"section"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"julia> include(\"PALEO_examples_sediment_Sisotopes.jl\")","category":"page"},{"location":"collated_examples/boudreau1996/README/","page":"Sediment Boudreau (1996) test cases","title":"Sediment Boudreau (1996) test cases","text":"three sediment column example as PALEO_examples_sediment.jl, with sulphur isotopes enabled and low (1 mM) oceanfloor [SO4] to illustrate Rayleigh fractionation within the sediment column as [SO4] becomes limiting.","category":"page"},{"location":"PALEOsediment_Reactions/#PALEOsediment-Reactions","page":"PALEOsediment Reactions","title":"PALEOsediment Reactions","text":"","category":"section"},{"location":"PALEOsediment_Reactions/#Sediment","page":"PALEOsediment Reactions","title":"Sediment","text":"","category":"section"},{"location":"PALEOsediment_Reactions/","page":"PALEOsediment Reactions","title":"PALEOsediment Reactions","text":"CurrentModule = PALEOsediment.Sediment","category":"page"},{"location":"PALEOsediment_Reactions/","page":"PALEOsediment Reactions","title":"PALEOsediment Reactions","text":"SedimentTransport.ReactionSedimentTransport\nSedimentBioRates.ReactionSedimentBioRates","category":"page"},{"location":"PALEOsediment_Reactions/#PALEOsediment.Sediment.SedimentTransport.ReactionSedimentTransport","page":"PALEOsediment Reactions","title":"PALEOsediment.Sediment.SedimentTransport.ReactionSedimentTransport","text":"SedimentTransport\n\nSediment transport for n x 1D sediment columns.\n\nA grid with n columns is created in a Domain sediment, bounded at the top by Domain oceanfloor and at the base by Domain sedimentfloor.  Boundary cells at the sediment surface are therefore in subdomain sediment.oceanfloor.\n\nThe number of columns, and column area, water depth, accumulation rate, and porosity are set by oceanfloor Variables. Bioturbation and bioirrigation rates should be supplied on the sediment grid eg by PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates.\n\nEach component <totalname> to be transported should be defined by a source-minus-sink flux <totalname>_sms, and one or more concentration Variables with names of form <rootnameN>_conc. Solute concentration Variables are identified by attributes vphase == VP_Solute and advect == true, and are transported by diffusion, bioturbation and bioirrigation, and advection. If the totalname attribute is present, then this is used to define the appropriate <totalname>_sms flux (allowing multiple species concentrations with different transport properties or phases for a single <totalname>), otherwise <totalname> is assumed to be the same as <rootnameN>. Transport fluxes are then accumulated into <totalname>_sms. Species-specific solute diffusivities are calculated based on the attribute diffusivity_speciesname of the <rootnameN>_conc solute Variables, which should be one of the  names available from PALEOaqchem.MolecularDiffusion.create_solute_diffusivity_func Similarly, solid phase concentration Variables are transported by bioturbation and advection, and identified by attribute vphase == VP_Solid and advect == true.\n\nOceanfloor solute fluxes should be defined in the Domain fluxOceanfloor, with names fluxOceanfloor.soluteflux_<totalname>. Input particulate fluxes should be added by the fluxOceanfloor flux coupler to the surface sediment cells by linking to Variables in the sediment.oceanfloor subdomain, ie to sediment.oceanfloor.<totalname>_sms.\n\nBurial fluxes are the base of the sediment columns are defined in the Domain fluxOceanBurial.\n\nParameters\n\nL[Float64]=0.15  (m), default_value=0.15, description=\"depth of sediment column\"\nncellspercol[Int64]=60, default_value=60, description=\"number of cells per column\"\nf_grid[String]=\"linear\", default_value=\"linear\", allowed_values=[\"linear\", \"quadratic\"], description=\"vertical grid transformation\"\ngrid_eta[Float64]=NaN  (m), default_value=NaN, description=\"length scale for vertical grid transformation\"\nf_porosity[String]=\"Const\", default_value=\"Const\", allowed_values=[\"Const\", \"ExpAtten\"], description=\"functional form for porosity vs depth\"\nzpor[Float64]=0.1  (m), default_value=0.1, description=\"lengthscale for porosity if f_porosity=ExpAtten\"\nzdbl[Float64]=0.0004  (m), default_value=0.0004, description=\"diffusive boundary layer thickness at sediment - water column interface\"\nw_solute[Bool]=false, default_value=false, description=\"true to assume w_solute = w_solid at base of column, false to set w_solute=0.0 (ie zero solute velocity at great depth)\"\n\nMethods and Variables for default Parameters\n\ndo_sediment_phys\nvolume  (m^3), VT_ReactDependency, description=\"volume of sediment cells\"\nvolume_total  (m^3), VT_ReactDependency, description=\"total volume of sediment cells\"\nAbox  (m^2), VT_ReactDependency, description=\"horizontal area of box\"\nzupper  (m), VT_ReactDependency, description=\"depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m\"\nzlower  (m), VT_ReactDependency, description=\"depth of lower surface of box (m)\"\nzmid  (m), VT_ReactDependency, description=\"mean depth of box\"\npressure  (dbar), VT_ReactDependency, description=\"sediment pressure\"\nrho_ref  (kg m^-3), VT_ReactDependency, description=\"density conversion factor\"\noceanfloor_Afloor –> oceanfloor.Afloor  (m^2), VT_ReactDependency, description=\"horizontal area of seafloor at sediment surface\"\noceanfloor_zfloor –> oceanfloor.zfloor  (m), VT_ReactDependency, description=\"depth of ocean floor (m, -ve)\"\noceanfloor_temp –> oceanfloor.temp  (K), VT_ReactDependency, description=\"oceanfloor temperature\"\noceanfloor_sal –> oceanfloor.sal  (psu), VT_ReactDependency, description=\"oceanfloor salinity\"\noceanfloor_phi –> oceanfloor.phi  (), VT_ReactDependency, description=\"sediment surface porosity\"\n[oceanfloor_phimin] –> oceanfloor.phimin  (), VT_ReactDependency, description=\"sediment porosity at infinite depth\"\noceanfloor_w_accum –> oceanfloor.w_accum  (m yr-1), VT_ReactDependency, description=\"sediment accumulation rate (+ve)\"\nphi  (), VT_ReactProperty, description=\"porosity (volume fraction of solute phase)\"\nphi_solid  (), VT_ReactProperty, description=\"1.0 - porosity (volume fraction of solid phase)\"\nvolume_solute  (m^3), VT_ReactProperty, description=\"solute volume of sediment cells\"\nvolume_solid  (m^3), VT_ReactProperty, description=\"solid volume of sediment cells\"\ntemp  (Kelvin), VT_ReactProperty, description=\"sediment temperature\"\nsal  (psu), VT_ReactProperty, description=\"sediment salinity\"\nw_solid  (m yr-1), VT_ReactProperty, description=\"solid phase advection velocity (downwards is -ve)\"\nw_solute  (m yr-1), VT_ReactProperty, description=\"solute phase advection velocity (downwards is -ve)\"\nDfac  (), VT_ReactProperty, description=\"tortuoisity-dependent multiplier for solute diffusivity\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOsediment_Reactions/#PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates","page":"PALEOsediment Reactions","title":"PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates","text":"ReactionSedimentBioRates\n\nCalculate sediment bioturbation and bioirrigation rates (for use by sediment transport).\n\nLiterature compilation of functional forms for dependency on oceanfloor oxygen, Corg flux, and depth within sediment, from (Boudreau, 1996), (Archer, 2002), (Arndt et al., 2011), (Dale et al., 2015)\n\nParameters\n\nf_bioTurbRate[String]=\"Prescribed\", default_value=\"Prescribed\", allowed_values=[\"Prescribed\", \"CorgArcher2002\"], description=\"functional form for bioturbation max rate\"\nf_bioTurbDepth[String]=\"ConstCutoff\", default_value=\"ConstCutoff\", allowed_values=[\"ConstCutoff\", \"ErfcCutoff\", \"ExpCutoff\"], description=\"functional form of bioturbation rate with depth in sediment\"\nf_bioIrrigDepth[String]=\"ConstCutoff\", default_value=\"ConstCutoff\", allowed_values=[\"ConstCutoff\", \"ErfcCutoff\", \"ExpCutoff\"], description=\"functional form of bioirrigation rate with depth in sediment\"\nf_bioO2[String]=\"None\", default_value=\"None\", allowed_values=[\"None\", \"MM\", \"Dale2015\"], description=\"functional form of bioturbation and bioirrigation sensitivity to oceanfloor oxygen\"\nbioO2halfmax[Float64]=0.02  (mol m-3), default_value=0.02, description=\"oceanfloor [O2] for 50% decrease in bioturbation/bioirrigation\"\nbioO2decreaserate[Float64]=0.012  (mol m-3), default_value=0.012, description=\"oceanfloor [O2] sharpness of decrease in bioturbation/bioirrigation\"\n\nMethods and Variables for default Parameters\n\nMETHODS PALEOboxes.DocStrings.Methods(:methods_do) exception: ErrorException(\"PALEOsediment.Sediment.SedimentBioRates.ReactionSedimentBioRates register_methods! not implemented\")\n\n\n\n\n\n","category":"type"},{"location":"References/#References","page":"References","title":"References","text":"","category":"section"},{"location":"References/","page":"References","title":"References","text":"Archer, D. E. (2002). A model of suboxic sedimentary diagenesis suitable for automatic tuning and gridded global domains. Global Biogeochemical Cycles 16.\n\n\n\nArndt, S.; Regnier, P.; Goddéris, Y. and Donnadieu, Y. (2011). GEOCLIM reloaded (v 1.0): a new coupled earth system model for past climate change. Geoscientific Model Development 4, 451–481.\n\n\n\nBoudreau, B. P. (1996). A method-of-lines code for carbon and nutrient diagenesis in aquatic sediments. Computers & Geosciences 22, 479–496.\n\n\n\nDale, A. W.; Nickelsen, L.; Scholz, F.; Hensen, C.; Oschlies, A. and Wallmann, K. (2015). A revised global estimate of dissolved iron fluxes from marine sediments. Global Biogeochemical Cycles 29, 691–707.\n\n\n\nVan Cappellen, P. and Wang, Y. (1996). Cycling of iron and manganese in surface sediments; a general theory for the coupled transport and reaction of carbon, oxygen, nitrogen, sulfur, iron, and manganese. American Journal of Science 296, 197–243.\n\n\n\n","category":"page"},{"location":"indexpage/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"indexpage/","page":"Index","title":"Index","text":"","category":"page"},{"location":"collated_examples/transport_tests/README/#Sediment-transport-test-cases","page":"Sediment transport test cases","title":"Sediment transport test cases","text":"","category":"section"},{"location":"collated_examples/transport_tests/README/","page":"Sediment transport test cases","title":"Sediment transport test cases","text":"include(\"PALEO_transport_sediment.jl\")","category":"page"},{"location":"collated_examples/transport_tests/README/","page":"Sediment transport test cases","title":"Sediment transport test cases","text":"Test solid and solute transport in sediment column with varying porosity.","category":"page"},{"location":"Equations_Units/#Equations,-units,-and-PALEO-variable-name-conventions","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"Conservation equations for a species C in PALEO are defined in terms of moles of C per cell, where a cell volume (m^3) is partitioned into solute volume volume_solute = volume * phi (m^3) and a solid volume volume_solid = volume * phi_solid (m^3) where phi is the porosity and phi_solid = 1 - phi.","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"The conservation equation for a species C is then:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"dC/dt = C_sms   (mol yr-1)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"A solute concentration is defined as moles per m^3 of solute:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"C_conc = C / volume_solute  (mol m-3)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"and a solid concentration as moles per m^3 of solid:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"C_conc = C /  volume_solid (mol m-3)","category":"page"},{"location":"Equations_Units/#Reservoir-configuration","page":"Equations, units, and PALEO variable name conventions","title":"Reservoir configuration","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"For the usual case where C is defined by a ReactionReservoir or ReactionReservoirTotal, the yaml file configuration for a solute reservoir should include:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"                variable_links:\n                    volume:                         volume_solute   # relink volume so C_conc is solute concentration\n                variable_attributes: \n                    R_conc:vphase:                  VP_Solute      # label phase for transport etc","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"and the yaml file for a solid reservoir should include:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"                variable_links:\n                    volume:                         volume_solid   # relink volume so C_conc is solid-phase concentration\n                variable_attributes: \n                    R_conc:vphase:                  VP_Solid      # label phase for transport etc","category":"page"},{"location":"Equations_Units/#Rate-laws-for-reactions","page":"Equations, units, and PALEO variable name conventions","title":"Rate laws for reactions","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"PALEO units for rate constants for both solute and solid-phase first-order rate constants are (mol m-3)-1 yr-1.  This follows the convention in (Van Cappellen and Wang, 1996) (see Table 3 and equations 72 and 76).","category":"page"},{"location":"Equations_Units/#Reactions-between-solute-species","page":"Equations, units, and PALEO variable name conventions","title":"Reactions between solute species","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"For a hypothetical reaction between solutes A, B, C","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"A + B -> C","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"with first-order kinetics hence reaction rate","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"R_abc = k_abc * A_conc * B_conc   (mol (m^3 solute)-1 yr-1)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"k_abc has units (mol m-3)-1 yr-1, where concentration units are mol (m^3 solute)-1.  A_sms and similar contributions to tracer conservation equations are then:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"A_sms = - R_abc * volume_solute  (mol yr-1)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"and similar for B_sms, C_sms.","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"This means that PALEO Reactions implementing eg secondary redox between solute species usually will require the yaml configuration to include:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"                variable_links:\n                    volume:                         volume_solute   # relink volume so volume refers to solute volume","category":"page"},{"location":"Equations_Units/#Reactions-between-a-solid-phase-species-and-solute-phase-species","page":"Equations, units, and PALEO variable name conventions","title":"Reactions between a solid phase species and solute phase species","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"For a hypothetical reaction between solid-phase species D and solute B","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"D + B -> C","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"with first-order kinetics hence reaction rate","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"R_dbc = k_dbc * D_conc * B_conc   (mol (m^3 solid)-1 yr-1)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"k_dbc (still) has units (mol m-3)-1 yr-1, where concentration units for D are mol (m^3 solid)-1 and for B are mol (m^3 solute)-1.  D_sms and similar contributions to tracer conservation equations are then:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"D_sms = - R_dbc * volume_solid  (mol yr-1)","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"and similar for B_sms, C_sms.","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"This means that PALEO Reactions implementing eg secondary redox between a solid and solute species usually will require the yaml configuration to include:","category":"page"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"                variable_links:\n                    volume:                         volume_solid   # relink volume so volume refers to solid phase volume","category":"page"},{"location":"Equations_Units/#Unit-conversions-for-rate-laws","page":"Equations, units, and PALEO variable name conventions","title":"Unit conversions for rate laws","text":"","category":"section"},{"location":"Equations_Units/","page":"Equations, units, and PALEO variable name conventions","title":"Equations, units, and PALEO variable name conventions","text":"[mol/l]-1 [yr-1]    = [mol m-3] [yr-1]\n1 M-1 yr-1          = 1e-3  (mol m-3)-1 yr-1\n1 M-0.5 yr-1        = 3.16e-2 (mol m-3)^-0.5 yr-1\n\n1 M-1 h-1           = 8.766  (mol m-3)-1 yr-1\n1 M-0.5 h-1         = 277.2  (mol m-3)^-0.5 yr-1","category":"page"},{"location":"#PALEOsediment.jl-documentation","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"","category":"section"},{"location":"#Installation-and-running-the-examples","page":"PALEOsediment.jl documentation","title":"Installation and running the examples","text":"","category":"section"},{"location":"#Installation","page":"PALEOsediment.jl documentation","title":"Installation","text":"","category":"section"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"NB: requires Julia 1.9 or later.  To check the Julia version:","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"julia> versioninfo()","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"Clone this github repository to local directory PALEOsediment: from a linux bash prompt or a Windows terminal,","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"$ git clone https://github.com/PALEOtoolkit/PALEOsediment.jl.git PALEOsediment","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"Start julia and navigate to the PALEOsediment/examples folder, and run setup.jl to configure the PALEOsediment/examples Julia environment to use the local (downloaded) version of the PALEOsediment package:","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"julia> cd(\"PALEOsediment/examples\")\njulia> include(\"setup.jl\") # use the local version of PALEOsediment packages to allow local modifications","category":"page"},{"location":"#Running-a-minimal-example","page":"PALEOsediment.jl documentation","title":"Running a minimal example","text":"","category":"section"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"Start julia and navigate to the PALEOsediment folder, then:","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"julia> cd(\"examples/boudreau1996\")\njulia> import Pkg\njulia> Pkg.activate(\"..\") # use the PALEOsediment/examples environment\n\njulia> include(\"PALEO_examples_sediment_x10.jl\")","category":"page"},{"location":"#Using-PALEOsediment-Reactions-from-other-models","page":"PALEOsediment.jl documentation","title":"Using PALEOsediment Reactions from other models","text":"","category":"section"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"The PALEO Reactions comprising the PALEOsediment models are available when the registered PALEOsediment package is loaded (without downloading the repository), ie","category":"page"},{"location":"","page":"PALEOsediment.jl documentation","title":"PALEOsediment.jl documentation","text":"julia> Pkg.add(\"PALEOsediment\")\njulia> import PALEOsediment","category":"page"}]
}
