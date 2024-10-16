isempty(ARGS) && throw(AssertionError("Please provide the path to your config.jl."))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your config.jl."))

config_file = abspath(ARGS[1])

isfile(config_file) || throw(AssertionError("Cannot find a file at $config_file."))

include(config_file)

project_path = dirname(config_file)

isfile(annotation_file) || (annotation_file = joinpath(project_path, annotation_file))
isfile(annotation_file) || throw(AssertionError("Cannot find $annotation_file. Please edit config.jl."))
isfile(genome_file) || (genome_file = joinpath(project_path, genome_file))
isfile(genome_file) || throw(AssertionError("Cannot find $genome_file. Please edit config.jl."))

isdir(data_folder) || (data_folder = joinpath(project_path, data_folder))
isdir(data_folder) || throw(AssertionError("Cannot find a directory called $data_folder. Please edit config.jl."))

using Pkg

Pkg.activate(joinpath(@__DIR__, "ChimericAnalysis"))
Pkg.instantiate()

using ChimericAnalysis, RNASeqTools

fastqs = (is_paired_end & !is_interleaved_paired_end) ?
    PairedSingleTypeFiles([(joinpath(data_folder, sname*suffix_read1*file_type), joinpath(data_folder, sname*suffix_read2*file_type)) for (sname, _) in samplename_condition]) :
    SingleTypeFiles([joinpath(data_folder, sname*file_type) for (sname, _) in samplename_condition])
filesexist(fastqs)

genome = Genome(genome_file)

types = unique(Vector{String}(vcat([srna_type, cds_type, rrna_type, trna_type], filter_types, additional_types)))
features = Features(annotation_file, types; name_keys=name_keys)
if autocomplete_utrs
    add5utrs!(features; cds_type=cds_type, utr_type=fiveutr_type, utr_length=autocomplete_utr_length)
    add3utrs!(features; cds_type=cds_type, utr_type=threeutr_type, utr_length=autocomplete_utr_length)
else
    merge!(features, Features(annotation_file, [threeutr_type, fiveutr_type]; name_keys=name_keys))
end

if autocomplete_igrs
    addigrs!(features; igr_type=igr_type)
else
    merge!(features, Features(annotation_file, [igr_type]; name_keys=name_keys))
end

merge_utrs_and_cds && (features = mergetypes(features, cds_type, threeutr_type, fiveutr_type, merge_type))

bp_parameters = (AU_score, GC_score, GU_score, bp_mismatch_penalty, bp_gap_open_penalty, bp_gap_extend_penalty)

trimmed = skip_trimming ? fastqs : trim_fastp(fastqs; fastp_bin=fastp_bin, min_length=min_length, average_window_quality=average_window_quality, deduplicate=deduplicate)

bams = align_mem(trimmed, genome;
    bwamem_bin=bwa_mem2_bin,
    samtools_bin=samtools_bin,
    is_interleaved_paired_end=is_interleaved_paired_end,
    min_score=min_alignment_score,
    match=match_score,
    mismatch=mismatch_penalty,
    gap_open=gap_open_penalty,
    gap_extend=gap_extend_penalty,
    clipping_penalty=clipping_penalty,
    unpair_penalty=unpair_penalty,
    unpair_rescue=unpair_rescue,
    min_seed_len=min_seed_len,
    reseeding_factor=reseeding_factor,
    sort_bam=sort_and_index_bam,
    threads=threads)

results_path = mkpath(joinpath(project_path, "results"))
conditions = Dict(c => [i for (i, info) in enumerate(samplename_condition) if info[2]===c] for c in unique([t[2] for t in samplename_condition]))
chimeric_analysis(features, bams, results_path, conditions, genome;
    filter_types=filter_types,
    min_distance=min_distance,
    prioritize_type=prioritize_srna ? srna_type : nothing,
    overwrite_type=igr_type,
    is_reverse_complement=is_reverse_complement,
    is_paired_end=is_paired_end,
    min_reads=min_reads,
    max_fisher_fdr=max_fisher_fdr,
    max_bp_fdr=max_bp_fdr,
    max_ligation_distance=max_ligation_distance,
    check_interaction_distances=bp_interval,
    include_read_identity=include_orientation,
    include_singles=include_singles,
    allow_self_chimeras=allow_self_chimeras,
    min_self_chimera_distance=min_self_chimera_distance,
    bp_parameters=bp_parameters,
    n_genome_samples=n_genome_samples,
    shift_weight=bp_shift_weight,
    keep_ints_without_ligation=keep_ints_without_ligation,
    filter_name_queries=filter_name_queries,
    join_pvalues_method=Symbol(combine_pvalues_method),
    min_mapping_quality=min_mapping_quality)
