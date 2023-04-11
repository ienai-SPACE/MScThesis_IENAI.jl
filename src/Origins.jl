"""   function: Origins

GOAL: 
     - Define the coordinates of the origins of the rays
     - The origins are placed on a plane, perpendicular to the velocity vector
     - The source of origins has a circular shape

INPUTS:
     - Vrel
     - distance: distance from body center to perpendicular plane from where the rays iniciate (the plane should be outside of the satellite body)  
     - rmax: radius of the ray source on the perpendicular plane 
OUTPUTS:
      - O         :: coordinates on the perpendicular plane from where the rays iniciate
      - Norig     :: number of ray origins

      PENDING:
      *unify the 'Areas' module
      *verify against the spherical tri-mesh
"""

function RayOrigins(dir, rmax, distance)


    Dangle= 20;           #Delta angle for the polar definition of origins

    #------pre-allocation-------------------
    xt = zeros(rmax,Int(360/Dangle));
    yt = zeros(rmax,Int(360/Dangle));
    zt = zeros(rmax,Int(360/Dangle));
    Xg = zeros(rmax*Int(360/Dangle),3);
    #---------------------------------------
        
    λ = atan(dir[2]/dir[1]); θ = asin(dir[3]/norm(dir));
    

    #rotation matrix: from the local ref. frame at the center of the satellite to the center of the source ref. frame
    R = [-sin(λ) cos(λ) 0;-sin(θ)*cos(λ) -sin(θ)*sin(λ) cos(θ);cos(θ)*cos(λ) cos(θ)*sin(λ) sin(θ)];
    
    x0_v = dir*distance;
    
    #create origin coordinates on the plane's local ref. frame
    for gg ∈ 0:Int((360/Dangle-1))  
        for rr ∈ 1:rmax
        xt[rr,gg+1] = rr*cos(gg*Dangle*π/180);
        yt[rr,gg+1] = rr*sin(gg*Dangle*π/180);
        zt[rr,gg+1]=0;
        end
    end
    
    #transform origin coordinates into the local satellite ref. frame    
    counter = 0;
    for gg ∈ 0:Int((360/Dangle -1))
        for rr ∈ 1:rmax
            counter += 1;
            Xg[counter,:] = transpose(R)*[xt[rr,gg+1];yt[rr,gg+1];zt[rr,gg+1]]+x0_v;
        end
    end
    
    Norig = rmax*Int(360/Dangle);
    return Xg, Norig
end