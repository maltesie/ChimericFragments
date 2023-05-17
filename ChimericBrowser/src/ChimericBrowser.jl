module ChimericBrowser

using PrecompileTools

using Dash, NetworkLayout, LightGraphs, GeometryBasics, Packing, DataStructures, Random, MultipleTesting
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

#@setup_workload begin
#
#    config_file = joinpath(assets_path, "precompile", "config.jl")
#    include(config_file)
#    project_path = dirname(config_file)
#    types = Dict(
#        "srna"=>srna_type,
#        "cds"=>merge_utrs_and_cds ? merge_type : cds_type,
#        "5utr"=>fiveutr_type,
#        "3utr"=>threeutr_type
#    )
#    param_dict::Vector{Pair{String,String}} = ["datasets" => "test",]
#    bp_parameters = (AU_score, GC_score, GU_score, bp_mismatch_penalty, bp_gap_open_penalty, bp_gap_extend_penalty)
#    check_interaction_distances = bp_interval
#    results_folder = joinpath(project_path, "results")
#    genome_file = joinpath(project_path, genome_file)
#    shift_weight = bp_shift_weight
#
#    @compile_workload begin
#
#        interactions, genome_info, genome = load_data(results_folder, genome_file, min_reads, max_fisher_fdr, max_bp_fdr)
#
#        scores = Dict(
#            (DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
#            (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
#            (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3]
#        )
#        model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
#                                                                gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])
#
#        seq_length = check_interaction_distances[1]-check_interaction_distances[2]
#        genomeseq = foldl(*, values(genome))
#        cgenomeseq = BioSequences.complement(genomeseq)
#        rgenomeseq = BioSequences.reverse(genomeseq)
#        rcgenomeseq = BioSequences.complement(rgenomeseq)
#        genome_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
#            i1 % 2 == 0 ? view(genomeseq, i1:i1+seq_length-1) : view(rcgenomeseq, i1:i1+seq_length-1),
#            i2 % 2 == 0 ? view(rgenomeseq, i2:i2+seq_length-1) : view(cgenomeseq, i2:i2+seq_length-1),
#            model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genomeseq)-seq_length), (n_genome_samples,2)))]
#        )
#
#        randseq = BioSequences.randdnaseq(length(genomeseq))
#        randseq_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
#            view(randseq, i1:i1+seq_length-1),
#            view(randseq, i2:i2+seq_length-1),
#           model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genomeseq)-seq_length), (n_genome_samples,2)))]
#        )
#
#        app = dash(assets_folder=assets_path)
#        app.layout = browser_layout(sort([k for k in keys(interactions)]), genome_info, stylesheet(types), min_reads, max_fisher_fdr, max_bp_fdr)
#
#        update_selection_callback!(app, interactions, seq_length, param_dict)
#        update_dataset_callback!(app, interactions, min_reads)
#        update_plots_callback!(app, interactions, randseq_model_ecdf, genome_model_ecdf)
#        update_selected_element_callback!(app, genome, interactions, check_interaction_distances, model)
#        click_cyto_button_callback!(app)
#        click_table_button_callback!(app, interactions)
#        click_add_node_callback!(app)
#    end
#end

end