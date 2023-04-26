module ChimericBrowser

using Dash, NetworkLayout, LightGraphs, GeometryBasics, Packing, DataStructures, Random
using BioAlignments, BioSymbols, FASTX, CSV, JLD2, DataFrames, StatsBase, PlotlyJS
import BioSequences

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))
const assets_path = realpath(joinpath( @__DIR__, "..", "assets"))
const rng = MersenneTwister(1234)

include("deps.jl")
include("data.jl")
include("cytostyle.jl")
include("plots.jl")
include("layout.jl")
include("callbacks.jl")
include("browser.jl")

end