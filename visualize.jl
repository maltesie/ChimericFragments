isempty(ARGS) && throw(AssertionError("Please provide the path to your config.jl."))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your config.jl."))

config_file = abspath(ARGS[1])

isfile(config_file) || throw(AssertionError("Cannot find a file at $config_file."))

include(config_file)

project_path = dirname(config_file)

isdir(joinpath(project_path, "results")) || throw(AssertionError("Cannot find results folder in the specified project folder. Please run analyze.jl first."))

isfile(joinpath(project_path, genome_file)) || throw(AssertionError("Cannot find a valid file with the filename $genome_file. Please edit config.jl."))

import Pkg

Pkg.activate("ChimericBrowser")
Pkg.instantiate()

using ChimericBrowser

types = Dict(
    "srna"=>srna_type,
    "cds"=>merge_utrs_and_cds ? merge_type : cds_type,
    "5utr"=>fiveutr_type,
    "3utr"=>threeutr_type
)

param_dict::Vector{Pair{String,String}} = [
    "datasets" => join(unique(v[2] for v in samplename_condition), ", "),
    "min. distance for chimeric classification:" => "$min_distance",
    "max. ligation distance:"=>"$max_ligation_distance",
    "min. reads per interaction:" => "$min_reads",
    "max. basepairing fdr:"=>"$max_bp_fdr",
    "max. Fisher's exact fdr:"=>"$max_fisher_fdr",
    "use order on read for Fisher's exact test:" => include_orientation ? "yes" : "no",
    "use single count for Fisher's exact test:" => include_singles ? "yes" : "no",
    "Fisher's test tail:"=>fisher_exact_tail,
    "self-chimeras included:" => allow_self_chimeras ? "yes" : "no",
    "autocompleted UTRs" => autocomplete_utrs ? "yes, with $autocomplete_utr_length nt max. length" : "no",
    "merged UTRs into CDS:" => merge_utrs_and_cds ? "yes" : "no",
]

bp_parameters = (AU_score, GC_score, GU_score, bp_mismatch_penalty, bp_gap_open_penalty, bp_gap_extend_penalty)

@info "Initializing visualization for datasets $(join(unique(v[2] for v in samplename_condition), ", "))."

chimeric_browser(joinpath(project_path, "results"), joinpath(project_path, genome_file), types, min_reads, max_fisher_fdr, max_bp_fdr,
    address, port, bp_interval, param_dict, bp_parameters, n_genome_samples, bp_shift_weight)
