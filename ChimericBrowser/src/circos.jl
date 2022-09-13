function circos(; kwargs...)
    available_props = Symbol[:id, :config, :enableDownloadSVG, :enableZoomPan, :eventDatum, :layout, :loading_state, :selectEvent, :size, :style, :tracks]
    wild_props = Symbol[]
    return Component("circos", "Circos", "dash_bio", available_props, wild_props; kwargs...)
end