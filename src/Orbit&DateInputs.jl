"""
    OrbitandDate()

Define orbital elements and date to obtain relative velocity of the satellite and atmospheric SurfaceProps
    
# Output:
- `JD`          : Julian date
- `alt`         : altitude [m]
- `g_lat`       : [rad]
- `g_long`      : [rad]
- `f107A`
- `f107`
- `ap`
- `Vrel_v`      : relative velocity vector [m/s]
"""
function OrbitandDate()

    #Orbital elements
    a = 400 + R_E     #[km]
    e = 0
    #r_p = a * (1 - e);        #[km]
    #r_a = a * (1 + e);        #[km]
    h = a - R_E          #[km]   mean altitude of the orbit

    Vatm = 0.1        #[km/s] (h = 200-300 km) velocity of rotation of the atmosphere --> WEST to EAST
    Vsc = sqrt(2 * mu / (h + R_E) - mu / (a))          #[km/s] velocity of the spacecraft
    Vrel = (Vsc - Vatm) * 1000  #[m/s] relative velocity between spacecraft and atmosphere



    Vrel_v = [1.0, 0.0, 0.0] / norm([0.0, 0.0, 1.0]) * Vrel #this vectorial velocity depends on the attitude of the S/C, since Vrel_v is given wrt body ref.frame




    #-----------NRLMSISE-00 INPUTS: Julian date, geodetic coordinates, and solar/magnetic indices-------------------------------
    JD = date_to_jd(2000, 1, 4, 0, 0, 0)     #Julian Day [UTC].
    alt = h * 1000                                #Altitude [m].
    g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
    g_long = deg2rad(-45)                       #Geodetic longitude [rad].
    f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79                                   #Daily F10.7 flux for previous day.
    ap = 5.13                                   #Magnetic index.
    #----------------------------------------------------------------------------------------------------------------------------



    return JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v

end