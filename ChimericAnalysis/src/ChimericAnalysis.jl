module ChimericAnalysis

using CSV, JLD2
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, StatsBase
using GenomicFeatures, IterTools, DataFrames, LoggingExtras, RNASeqTools, BioSequences, BioAlignments

import Plots: plot, plot!, histogram, histogram!, title!, vline!, xlabel!, ylabel!, heatmap, @layout, grid
using Measures

export chimeric_analysis, mergetypes

include("chimeric.jl")
include("plots.jl")

end
