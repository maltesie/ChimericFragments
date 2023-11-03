# deps.jl defines function for loading external JavaScript components from the deps folder

# load the circos plot component from DashBio and make properties available
function circos(; kwargs...)
    available_props = Symbol[:id, :config, :enableDownloadSVG, :enableZoomPan, :eventDatum, :layout, :loading_state, :selectEvent, :size, :style, :tracks]
    wild_props = Symbol[]
    return Component("circos", "Circos", "dash_bio", available_props, wild_props; kwargs...)
end

# load the cytoscape component for drawing graphs from DashBio and make properties available
function cytoscape(; kwargs...)
    available_props = Symbol[:id, :autoRefreshLayout, :autolock, :autoungrabify, :autounselectify, :boxSelectionEnabled, :className,
                                :elements, :generateImage, :imageData, :layout, :maxZoom, :minZoom, :mouseoverEdgeData, :mouseoverNodeData,
                                :pan, :panningEnabled, :responsive, :selectedEdgeData, :selectedNodeData, :style, :stylesheet, :tapEdge,
                                :tapEdgeData, :tapNode, :tapNodeData, :userPanningEnabled, :userZoomingEnabled, :zoom, :zoomingEnabled]
    wild_props = Symbol[]
    return Component("cytoscape", "Cytoscape", "dash_cytoscape", available_props, wild_props; kwargs...)
end

# register external components and let Dash handle the communication between the component and ChimericFragments
function __init__()
    DashBase.register_package(
        DashBase.ResourcePkg(
            "dash_cytoscape",
            resources_path,
            version = "0.3.0",
            [
                DashBase.Resource(
                    relative_package_path = "dash_cytoscape_extra.min.js",
                    external_url = "https://unpkg.com/dash-cytoscape@0.3.0/dash_cytoscape/dash_cytoscape_extra.min.js",
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
