module ChimericAnalysis

using XAM, FASTX, CodecZlib, GFF3, BigWig, DelimitedFiles, BGZFStreams, CSV, XLSX
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, StatsBase
using ElasticArrays, IterTools, LightGraphs, DataFrames

export chimeric_analysis, add5utrs!, add3utrs!, addigrs!, Genome, Features, merge!, align_mem
export SingleTypeFiles, PairedSingleTypeFiles, check_files_exist

include("types.jl")
include("preprocess.jl")
include("files.jl")
include("sequence.jl")
include("annotation.jl")
include("alignment.jl")
include("chimeric.jl")

end
