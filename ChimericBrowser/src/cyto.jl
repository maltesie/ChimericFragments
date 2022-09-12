using Pkg
Pkg.activate("/home/abc/Workspace/RNASeqAnalysis/")

using Dash, DashHtmlComponents, DashCoreComponents

DashBase.register_package(
        DashBase.ResourcePkg(
            "dash_cytoscape",
            resources_path,
            version = version,
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

app = dash()

app.layout = html_div([
    cyto_cytoscape(
        id="cytoscape-two-nodes",
        layout=Dict("name" =>  "preset"),
        style=Dict("width" =>  "100%", "height" =>  "400px"),
        elements=[
            Dict("data" =>  Dict("id" =>  "one", "label" =>  "Node 1"), "position" =>  ("x" =>  75, "y" =>  75)),
            Dict("data" =>  Dict("id" =>  "two", "label" =>  "Node 2"), "position" =>  ("x" =>  200, "y" =>  200)),
            Dict("data" =>  Dict("source" =>  "one", "target" =>  "two"))
        ]
    )
])

run_server(app, "0.0.0.0", debug=true)
