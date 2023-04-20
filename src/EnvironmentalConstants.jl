#=Environmental constants: 

- Gravitational parameter of the Earth
- Boltzmann constant
- Avogadro's number
- Earth radius
- Semi-major axis
- Eccentricity
- Orbit perigee and apogee
- Mean altitude of the orbit
- Atmospheric, spacecraft, and relative velocity
- Oxygen concentration
- Mean molecular mass
- Thermal speed of the gas
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

#=
#Orbital elements
const a = 200 + R_E;     #[km]
const e = 0;
const r_p = a * (1 - e);        #[km]
const r_a = a * (1 + e);        #[km]


const h = a - R_E;          #[km]   mean altitude of the orbit


const Vatm = 0.1;        #[km/s] (h = 200-300 km) velocity of rotation of the atmosphere
const Vsc = sqrt(2 * mu / (h + R_E) - mu / (a));          #[km/s] velocity of the spacecraft
const Vrel = (Vsc - Vatm) * 1000;  #[m/s] relative velocity between spacecraft and atmosphere


#From MATLAB MSISE00 model:
# C = [O H N O2 N2 He]
C = @SVector[0.757036209368285, 1.03217482656703e-05, 0.00664130845103655, 0.00627275013672124, 0.228876415580270, 0.000755673189780810];

const P0 = 2e-5;         #[Pa] for 300km (5×10−6;10−6 for 400km and 500km, respectively) computed with NRLMSISE-00
#this should be calculated witht the atmospheric model

const mmean = 1;         #[kg] mean molecular mass
#this should be calculated witht the atmospheric model

const Ta = 943.83;       #[K] ambient temperature
#this should be calculated witht the atmospheric model

#const s = Vsc*1000/sqrt(2*kb*Ta/mmean);   # [-] thermal speed of the gas
=#

export mu, kb, R, NA, R_E