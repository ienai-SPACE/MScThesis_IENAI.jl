using JSON
import FilePathsBase
using FilePathsBase: /

"""
    load_material_properties()

Load .json file linking the facet index with its corresponding surface material

#OUTPUT
- `mesh_facets::Vector`     : vector of `Dict{String, Any}` containing the field `"material_index"` (indexing starts in 0). The field `"material_index"` corresponds to the position of the material in the vector `materials`
- `materials::Vector`       : vector of `Dict{String, Any}` containing the fields `"name"` and `"atomic mass"`
"""

function load_material_properties(path)
    # pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    # Load the materials and mesh facets from "materials.json" file
    json_data = open(path) do file
        read(file, String)
    end


    input_material_data = JSON.parse(json_data)
    materials = input_material_data["materials"]
    mesh_facets = input_material_data["mesh_facets"]

    # # Process each facet
    # for facet in mesh_facets
    #     material_index = facet["material_index"] + 1  # Julia uses 1-based indexing
    #     material_properties = materials[material_index]

    #     # Now you can use the material_properties as needed for each facet
    #     println("Facet with material index ", material_index, " has properties: ", material_properties)
    # end

    return mesh_facets, materials

end

function load_material_properties()
    pkg_path = (FilePathsBase.@__FILEPATH__() |> parent |> parent)
    load_material_properties(pkg_path / "test" / "inputs_models_data" / "facetMaterials.json")
end