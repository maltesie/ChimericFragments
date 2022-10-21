isempty(ARGS) && throw(AssertionError("Please provide the path to your config.jl."))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your config.jl."))

config_file = ARGS[1]

isfile(config_file) || throw(AssertionError("Cannot find a file at $config_file."))

include(config_file)

project_path = dirname(config_file)

isdir(joinpath(project_path, "results")) || throw(AssertionError("Cannot find results folder in the specified project folder. Please run analyze.jl first."))

isfile(joinpath(project_path, genome_file)) || throw(AssertionError("Cannot find a valid file with the filename $genome_file. Please edit config.jl."))

using ChimericBrowser

chimeric_browser(joinpath(project_path, "results"), joinpath(project_path, genome_file), srna_type)
