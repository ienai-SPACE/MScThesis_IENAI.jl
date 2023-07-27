using JSON
import FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# Load the materials and mesh facets from "materials.json" file
json_data = open(pkg_path / "test" / "inputs_models_data" / "facetMaterials.json") do file
    read(file, String)
end


input_material_data = JSON.parse(json_data)
materials = input_material_data["materials"]
mesh_facets = input_material_data["mesh_facets"]

# Process each facet
for facet in mesh_facets
    material_index = facet["material_index"] + 1  # Julia uses 1-based indexing
    material_properties = materials[material_index]

    # Now you can use the material_properties as needed for each facet
    println("Facet with material index ", material_index, " has properties: ", material_properties)
end
