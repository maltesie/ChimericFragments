module ChimericBrowser

using Dash, DataFrames, Random, CSV, FASTX, JLD2, NetworkLayout, LightGraphs, GeometryBasics, Packing, DataStructures
using BioAlignments: pairalign, AffineGapScoreModel, LocalAlignment, count_matches, alignment
using BioSequences: reverse_complement, LongDNA

export chimeric_browser

const resources_path = realpath(joinpath( @__DIR__, "..", "deps"))
const assets_path = realpath(joinpath( @__DIR__, "..", "assets"))
const rng = MersenneTwister(1234)

include("deps.jl")
include("data.jl")
include("cytostyle.jl")
include("layout.jl")
include("callbacks.jl")
include("browser.jl")

end