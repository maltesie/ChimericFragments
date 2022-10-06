using Pkg

isempty(ARGS) && throw(AssertionError("Please provide the path to your config.jl."))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your config.jl."))

config_file = ARGS[1]

isfile(config_file) || throw(AssertionError("Cannot find a file at $config_file."))

include(config_file)

project_path = dirname(config_file)

isfile(joinpath(project_path, annotation_file)) || throw(AssertionError("Cannot find a file with the filename $annotation_file. Please edit config.jl."))
isfile(joinpath(project_path, genome_file)) || throw(AssertionError("Cannot find a file with the filename $genome_file. Please edit config.jl."))

data_path = joinpath(project_path, data_folder)
isdir(joinpath(project_path, data_folder)) || throw(AssertionError("Cannot find a directory called $data_folder. Please edit config.jl."))

Pkg.activate(joinpath(@__DIR__, "ChimericAnalysis"))
Pkg.instantiate()

using ChimericAnalysis, BioSequences

fastqs = is_paired_end ? PairedSingleTypeFiles([(joinpath(data_path, sname*suffix_read1*".fastq.gz"), joinpath(data_path, sname*suffix_read2*".fastq.gz"))
                                                for (sname, _) in samplename_condition]) :
                            SingleTypeFiles([joinpath(data_path, sname*".fastq.gz") for (sname, _) in samplename_condition])
check_files_exist(fastqs)

genome = Genome(joinpath(project_path, genome_file))

types = vcat([srna_type, cds_type, rrna_type, trna_type], additional_types)
features = Features(joinpath(project_path, annotation_file), types)
if autocomplete_utrs
    add5utrs!(features; cds_type=cds_type, utr_type=fiveutr_type, utr_length=autocomplete_utr_length)
    add3utrs!(features; cds_type=cds_type, utr_type=threeutr_type, utr_length=autocomplete_utr_length)
else
    merge!(features, Features(joinpath(project_path, annotation_file), [threeutr_type, fiveutr_type]))
end

if autocomplete_igrs
    addigrs!(features; igr_type=igr_type)
else
    merge!(features, Features(joinpath(project_path, annotation_file), [igr_type]))
end

features = mergetypes(features, [threeutr_type, fiveutr_type, cds_type], "UTRS_CDS")

bams = align_mem(fastqs, genome; bwa_bin=bwa_mem2_bin, sam_bin=samtools_bin)

results_path = mkpath(joinpath(project_path, "results"))
conditions = Dict(c => [i for (i, info) in enumerate(samplename_condition) if info[2]===c] for c in unique([t[2] for t in samplename_condition]))
chimeric_analysis(features, bams, results_path, conditions; filter_types=filter_types, min_distance=min_distance, prioritize_type=prioritize_srna ? srna_type : nothing,
    overwrite_type=igr_type, is_reverse_complement=is_reverse_complement, min_reads=min_reads, max_fdr=max_fdr,
    include_read_identity=include_read_identity, include_singles=include_singles)
