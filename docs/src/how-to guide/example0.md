# Example 0: inputs formats

## Geometry input .obj file example

As mentioned in the *Input/Output discussion*, it is essential that all triangular facets follow the same wind-up scheme. For this, an *.obj* file exported from *Blender* (after wind-up scheme unification) is a fast method to ensure input correctness.

```
# Blender 3.5.1
# www.blender.org
mtllib T_Sat.mtl
o T_Sat
v 299.997986 -0.003553 0.005614
v 299.997986 49.999699 0.002687
v 300.000000 0.000582 74.999802
v 299.997986 49.997505 75.004303
v 299.998993 -0.000882 150.001999
v 300.001007 49.997013 150.001999
v 299.997986 0.001358 225.001007
v 299.998993 49.999817 225.001999
v 300.000000 -0.000996 300.001007
v 299.998993 49.998722 300.002991
v 299.998993 0.000465 375.001007
v 300.000000 49.999226 375.001007

...
```


## Material input .json file example

In order to simplify this functionality for the user, a widely used and self-explanatory format is selected to include the properties: .json.

```
{
    "materials": [
        {
            "name": "aluminum",
            "atomic mass": 26.9815
        },
        {
            "name": "monocrystalline silicon",
            "atomic mass": 32.12
        },
        {
            "name": "silicon",
            "atomic mass": 28.0855
        }
    ],
    "mesh_facets": [
        {
            "material_index": 0
        },
        {
            "material_index": 0
        },
        {
            "material_index": 1
        },
        {
            "material_index": 1
        },
        {
            "material_index": 0
        },
        {
            "material_index": 0
        },

               ...

    ]
}
```