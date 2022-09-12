DashBase.register_package(
    DashBase.ResourcePkg(
        "dashbio_circos",
        resources_path,
        version = "0.7.0",
        [
            DashBase.Resource(
                relative_package_path = "async-circos.js",
                external_url = "https://unpkg.com/dash-bio@0.7.0/dash_bio/async-circos.js",
                dynamic = nothing,
                async = :true,
                type = :js
            )
        ]
    )
)

function dash_circos(; kwargs...)
    available_props = Symbol[:id, :config, :enableDownloadSVG, :enableZoomPan, :eventDatum, :layout, :loading_state, :selectEvent, :size, :style, :tracks]
    wild_props = Symbol[]
    return Component("dashbio_circos", "DashBio", "dash_circos", available_props, wild_props; kwargs...)
end