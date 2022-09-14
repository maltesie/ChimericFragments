isempty(ARGS) && throw(AssertionError("Please provide the path to your project folder containing the results of the chimeric fragments analysis."))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your project folder."))
isdir(ARGS[1]) || throw(AssertionError("Please provide a valid path."))

data_path = ARGS[1]

isdir(joinpath(data_path, "results")) || throw(AssertionError("Cannot find results folder in the specified project folder. Please run analyze.jl first."))
isfile(joinpath(data_path, "config.jl")) || throw(AssertionError("Please add a config.jl to the specified folder."))

include(joinpath(data_path, "config.jl"))
isfile(joinpath(data_path, genome_file)) || throw(AssertionError("Cannot find a valid file with the filename $genome_file. Please edit config.jl."))

using Pkg
Pkg.activate(joinpath(@__DIR__, "ChimericBrowser"))
Pkg.instantiate()

using ChimericBrowser

chimeric_browser(joinpath(data_path, "results"), genome_file)