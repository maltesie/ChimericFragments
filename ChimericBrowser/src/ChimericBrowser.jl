module ChimericBrowser

using Dash

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))

include("circos.jl")
include("cyto.jl")
include("annotation.jl")
include("data.jl")
include("layout.jl")

function chimeric_browser(project_path::String)
    app = dash()
    app.layout = browser_layout
    run_server(app, "0.0.0.0", debug=true)
end

end