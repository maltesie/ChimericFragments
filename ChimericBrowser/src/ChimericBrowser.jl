module ChimericBrowser

using Dash, DataFrames, CSV, FASTX

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))
const assets_path = realpath(joinpath( @__DIR__, "..", "assets"))

include("deps.jl")
include("data.jl")
include("layout.jl")
include("callbacks.jl")

function chimeric_browser(results_folder::String, genome_file::String)

    interactions_dfs, singles_dfs, stats_dfs, gene_names, genome_info = load_data(results_folder, genome_file)

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions_dfs)]), [Dict("label"=>"test", "value"=>"test")], genome_info)

    update_selection_callback!(app, interactions_dfs)

    run_server(app, "0.0.0.0", debug=true)
end

end