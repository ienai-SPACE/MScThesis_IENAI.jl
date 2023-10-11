using Meshes
import WGLMakie as Mke

geom = Box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

# mesh = simplexify(boundary(geom))

# hexagon = Hexagon((0.0, 0.0), (1.0, 0.0), (1.0, 1.0),
#     (0.75, 1.5), (0.25, 1.5), (0.0, 1.0))

# mesh = discretize(hexagon, FanTriangulation()) |> viz

sphere = Sphere((0.0, 0.0, 0.0), 1.0)
mesh = discretize(sphere, RegularDiscretization(10, 10))

sphereMesh = discretize(geom, RegularDiscretization(10, 10))
mesh = refine(sphereMesh, TriRefinement())

mesh3 = discretize(geom, BoundaryDiscretizationMethod)

a = Base.convert(Hexahedron, geom)

hexa = Hexahedron((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))

mesh3 = discretize(boundary(geom)) |> viz
refinedMesh = refine(mesh3, TriSubdivision())

aa = Multi(geom)