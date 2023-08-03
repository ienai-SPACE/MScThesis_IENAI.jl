#=Environmental constants: 

- Gravitational parameter of the Earth
- Boltzmann constant
- Gas constant [J K^-1 mol^-1]
- Avogadro's number
- Earth radius
=#

using StaticArrays

#Gravitational parameter of the Earth
const mu = 3.98600e5; #[km^2/s]
#Boltzmann constant
const kb = 1.380649e-23;    #[m^2 kg s^-2 K^-1]
const R = 8.31446261815324; # Gas constant [J K^-1 mol^-1]
#Avogadro's number
const NA = 6.0221408e23;
#Earth radius
const R_E = 6371;   #[km]

export mu, kb, R, NA, R_E