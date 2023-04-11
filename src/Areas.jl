

using LinearAlgebra
 
include("MTalgorithm.jl")
include("Origins.jl")


""" Areas:
GOAL:
    - Obtaining all the areas of the triangular mesh elements that have been intercepted by the rays (areas are not projected onto the velocity directions)
        - The rays are originated on a plane which is perpendicular to the velocity vector, and far from the satellite body 


INPUT:
    - rmax            :: radius of the circular plane from where rays originate
    - distance        :: distance at which the circular plane is located (it should be out from the satellite body)
    - Vrel            :: relative velocity vector --> used to calculate the perpendicular plane to the satellite's velocity direction

OUTPUT:
    - OutTriangles    :: matrix storing the triangle index, area, and the angle between the normal and the velocity direction
        - Matrix{Float64} -> size: (3, number of intercepted triangles)
    - Aproj           :: projection of the intercepted triangular areas onto the velocity direction 
    - Aref            :: sum of all intercepted triangular areas
"""


rmax = 1;                            #radius of the circular plane from where rays originate
distance = 10;                       #distance at which the circular plane is located (it should be out from the satellite body)
Vrel = [7000 7000 7000];

#direction vector
dir = [-0.6988494037821777, -0.137916655763437, -0.7018464981007773];
#dir = -(Vrel/norm(Vrel))';       
#triangle vetices coordinates
triangles = [4.394897596952136 -1.3063587207678875 4.012655067923802 3.3442061435039823 0.7676371202562053 4.251153740481179 5.445214679778381 1.8739750984535304 5.439657623886108;
1.3410577524314113 0.6227469916184274 1.5604711027511295 1.4074299081230686 -0.5929514915580713 2.0821525791882842 2.850415168046063 0.5144988358467968 1.1637316942088742];
#Number of triangles
Ntri = size(triangles,1);



O, Norig = RayOrigins(dir, rmax, distance);        #coordinates of ray origins, number of origins

#------pre-allocation-------------------
index = 0;                                         #counter indicating the number of triangles intercepted by the same ray
intercept_dummy = zeros(Ntri,4);                   #triangle index, distance from origin to intercept, area of the triangle, angle between velocity vector and triangle's normal
triIntercept = zeros(3,Norig);                     #triangle index, area of the triangle, angle between velocity vector and triangle's normal
#---------------------------------------

for jj ∈ 1:Norig       #iterate over the set of ray origins
    for ii ∈ 1:Ntri    #iterate over all triangles
        
        #definition of the triangle vertices
        V1 = [triangles[ii,1],triangles[ii,2],triangles[ii,3]];
        V2 = [triangles[ii,4],triangles[ii,5],triangles[ii,6]];
        V3 = [triangles[ii,7],triangles[ii,8],triangles[ii,9]];

        area, flag, facing, u, v, t, γ_dir = MTalgorithm(O[jj,:], dir, V1, V2, V3);

        if flag == 0

        elseif flag == 1   #the triangle is intercepted by the ray
          
            global index
            index += 1;

            intercept_dummy[index,:] = [ii,t,area, γ_dir];

            if index > 1 #if more than one triangle intercepted by the same ray
                if abs(intercept_dummy[index,2])< abs(intercept_dummy[index-1,2])   #select the most forward triangle to be intercepted
                
                    triIntercept[:,jj] = [intercept_dummy[index,1] intercept_dummy[index,3] intercept_dummy[index,4]];   #triangle index, area of the triangle, angle between velocity vector and triangle's normal

                end
            else
                triIntercept[:,jj] = [intercept_dummy[index,1] intercept_dummy[index,3] intercept_dummy[index,4]];           #triangle index, area of the triangle, angle between velocity vector and triangle's normal
            end         
        end
    end

    global index
    index = 0;
    fill!(intercept_dummy, 0);

end


#------pre-allocation-------------------
OutTriangles = zeros(3,Norig);            #matrix storing the triangle index, area, and the angle between the normal and the velocity direction
triangleCount = 0;
indexList = zeros(Norig,1);               #list to store the indices of the triangles that are intercepted
#---------------------------------------

#check the triangles that have been intercepted and store their properties
for jj ∈ 1:Norig
    if Int(triIntercept[1,jj]) != 0
        if !(triIntercept[1,jj] ∈ indexList)
            global triangleCount
            triangleCount += 1;
            OutTriangles[:,triangleCount] = triIntercept[:,jj];
            indexList[triangleCount] = triIntercept[1,jj];
        end
    end
end

#eliminate non-intercepted triangle entries and size down the output matrix
OutTriangles = filter(!iszero, OutTriangles);
OutTriangles = reshape(OutTriangles,(3,Int(length(OutTriangles)/3)))

Aref = sum(OutTriangles[2,:])      #sum of all intercepted triangular areas

#------pre-allocation-------------------
Aproj = zeros(size(OutTriangles,2),1);
#---------------------------------------


for ii ∈ 1:size(OutTriangles,2)
    Aproj[ii] = abs(OutTriangles[2,ii]*cos(OutTriangles[3,ii]));
end
Aproj = sum(Aproj)                 #sum of all intercepted triangular projected areas