module ChimericAnalysis

using CSV, XLSX, JLD2
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, StatsBase
using GenomicFeatures, IterTools, DataFrames, LoggingExtras, RNASeqTools

export chimeric_analysis, mergetypes

include("chimeric.jl")

end
