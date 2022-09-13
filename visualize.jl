#isempty(ARGS) && throw(AssertionError("Please provide the path to your project folder containing the results of the chimeric fragments analysis."))
#length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your project folder."))
#isdir(ARGS[1]) || throw(AssertionError("Please provide a valid path."))

#data_path = ARGS[1]

#isdir(joinpath(data_path, "results")) || throw(AssertionError("Cannot find results folder in the specified project folder. Please run analyze.jl first."))
using Pkg
Pkg.activate(joinpath(@__DIR__, "ChimericBrowser"))
Pkg.instantiate()

using ChimericBrowser

test_datapaths = ["/home/abc/super.csv"]
test_function_categories = [Dict("label"=>"hallo", "value"=>"hallo")]
test_genome_info = ["chrom1"=>1000000, "chrom2"=>4000000]

chimeric_browser(test_datapaths, test_function_categories, test_genome_info)