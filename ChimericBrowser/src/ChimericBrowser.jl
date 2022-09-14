module ChimericBrowser

using Dash, DataFrames, CSV
using FASTX.FASTA

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))
const assets_path = realpath(joinpath( @__DIR__, "..", "assets"))

include("deps.jl")
include("data.jl")
include("layout.jl")
include("callbacks.jl")

function chimeric_browser(dataset_paths::Vector{String}, function_categories::Vector{Dict{String,String}}, genome_info::Vector{Pair{String,Int}})
    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(dataset_paths, function_categories, genome_info)
    run_server(app, "0.0.0.0", debug=true)
end

end