"""
    Gas{T}
    ------------------------------------

    Gas atomic weights

 - m::T    # molecular mass [g/mol]


    srfAtomProps{T}
    ------------------------------------

    Surface material atomic masses

 - mass::T      #mass of an atom [g]
 - amass::T     #atomic mass [u] 


    SRMatProps{T} 
    ------------------------------------

    Solar radiation material properties

 - ads::T     #adsorptance
 - emitt::T    #emittance


   other constants:
   ------------------------------------
   eta  : generator efficienty
   Tw   : wall temperature
   s_cr : specular reflectivity component
   s_cr : diffuse reflectivity component
"""


using StaticArrays

#=
const η = 0;   #generator efficienty

const Tw = 300; #[K] temperature at the wall        

#= Solar radiation pressure coefficients (for calculation of force of the satellite due to
solar radiation according to [Luthcke et al., 1997]) =#

const s_cr = 0.15; #specular reflectivity component
const s_cd = 0.25; #diffuse reflectivity component
=#


struct Gas{T}
    m::T    # molecular mass [g/mol]
end

const O = Gas(15.9994);
const H = Gas(1.0079);
const N = Gas(14.0067);
const O2 = Gas(31.9988);
const N2 = Gas(28.0134);
const He = Gas(2.0);
const m = @SVector[O, H, N, O2, N2, He];


struct srfAtomProps{T}
    mass::T      #mass of an atom [g]
    amass::T     #atomic mass [g/mol]
end

const Aluminum = srfAtomProps(4.48e-23, 26.9815);
const Silicon = srfAtomProps(4.664e-23, 28.0855);


struct SRMatProps{T}

    #struct to store solar radiation material properties

    ads::T     #adsorptance
    emitt::T    #emittance

end

const Polished_beryllium = SRMatProps(0.44, 0.01);
const Goldized_kapton = SRMatProps(0.25, 0.02);
const Aluminum_tape = SRMatProps(0.21, 0.04);
const Solar_cells_Si = SRMatProps(0.75, 0.82);

export eta, Tw, s_cr, s_cd, O, H, N, O2, N2, He, m, Aluminum, Silicon, Polished_beryllium, Goldized_kapton, Aluminum_tape, Solar_cells_Si