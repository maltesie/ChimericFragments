module ChimericBrowser

using Dash

export chimeric_browser, cyto_browser, circos_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))

include("circos.jl")
include("cyto.jl")
include("annotation.jl")
include("data.jl")
include("layout.jl")

function chimeric_browser()
    app = dash()
    app.layout = browser_layout
    run_server(app, "0.0.0.0", debug=true)
end

function __init__()
    DashBase.register_package(
        DashBase.ResourcePkg(
            "dash_cytoscape",
            resources_path,
            version = "0.3.0",
            [
                DashBase.Resource(
                    relative_package_path = "dash_cytoscape.min.js",
                    external_url = "https://unpkg.com/dash-cytoscape@0.3.0/dash_cytoscape/dash_cytoscape.min.js",
                    dynamic = nothing,
                    async = nothing,
                    type = :js
                )
            ]
        )
    )
    DashBase.register_package(
        DashBase.ResourcePkg(
            "dash_bio",
            resources_path,
            version = "0.7.0",
            [
                DashBase.Resource(
                    relative_package_path = "async-circos.js",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/async-circos.js",
                    dynamic = nothing,
                    async = :true,
                    type = :js
                ),
                DashBase.Resource(
                    relative_package_path = "async-circos.js.map",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/async-circos.js.map",
                    dynamic = true,
                    async = nothing,
                    type = :js
                ),
                DashBase.Resource(
                    relative_package_path = "bundle.js",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/bundle.js",
                    dynamic = nothing,
                    async = nothing,
                    type = :js
                ),
                DashBase.Resource(
                    relative_package_path = "bundle.js.map",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/bundle.js.map",
                    dynamic = true,
                    async = nothing,
                    type = :js
                ),
                DashBase.Resource(
                    relative_package_path = "dash_bio-shared.js",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/dash_bio-shared.js",
                    dynamic = nothing,
                    async = :true,
                    type = :js
                ),
                DashBase.Resource(
                    relative_package_path = "dash_bio-shared.js.map",
                    external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/dash_bio-shared.js.map",
                    dynamic = true,
                    async = nothing,
                    type = :js
                )
            ]
        )
    )
end

end