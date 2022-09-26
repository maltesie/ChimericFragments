module ChimericBrowser

using Dash, DataFrames, CSV, FASTX

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))
const assets_path = realpath(joinpath( @__DIR__, "..", "assets"))

include("deps.jl")
include("data.jl")
include("layout.jl")
include("callbacks.jl")
include("browser.jl")

end