

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

function cyto_cytoscape(; kwargs...)
    available_props = Symbol[:id, :autoRefreshLayout, :autolock, :autoungrabify, :autounselectify, :boxSelectionEnabled, :className,
                                :elements, :generateImage, :imageData, :layout, :maxZoom, :minZoom, :mouseoverEdgeData, :mouseoverNodeData,
                                :pan, :panningEnabled, :responsive, :selectedEdgeData, :selectedNodeData, :style, :stylesheet, :tapEdge,
                                :tapEdgeData, :tapNode, :tapNodeData, :userPanningEnabled, :userZoomingEnabled, :zoom, :zoomingEnabled]
    wild_props = Symbol[]
    return Component("cyto_cytoscape", "Cytoscape", "dash_cytoscape", available_props, wild_props; kwargs...)
end