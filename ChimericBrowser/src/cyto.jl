function cytoscape(; kwargs...)
    available_props = Symbol[:id, :autoRefreshLayout, :autolock, :autoungrabify, :autounselectify, :boxSelectionEnabled, :className,
                                :elements, :generateImage, :imageData, :layout, :maxZoom, :minZoom, :mouseoverEdgeData, :mouseoverNodeData,
                                :pan, :panningEnabled, :responsive, :selectedEdgeData, :selectedNodeData, :style, :stylesheet, :tapEdge,
                                :tapEdgeData, :tapNode, :tapNodeData, :userPanningEnabled, :userZoomingEnabled, :zoom, :zoomingEnabled]
    wild_props = Symbol[]
    return Component("cytoscape", "Cytoscape", "dash_cytoscape", available_props, wild_props; kwargs...)
end