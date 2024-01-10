module ChimericAnalysis

using CSV, JLD2
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, StatsBase
using GenomicFeatures, IterTools, DataFrames, PrettyTables, LoggingExtras, RNASeqTools, BioSequences, BioAlignments

export chimeric_analysis, mergetypes

include("mergedreads.jl")
include("interactions.jl")
include("pipeline.jl")

end
