using Pkg

isempty(ARGS) && throw(AssertionError("Please provide the path to your RIL-Seq data folder containing the two fastq files and a config.jl"))
length(ARGS) > 1 && throw(AssertionError("Please provide only one path to your RIL-seq folder."))
isdir(ARGS[1]) || throw(AssertionError("Please provide a valid path."))

data_path = ARGS[1]

isfile(joinpath(data_path, "config.jl")) || throw(AssertionError("Please add a config.jl to the specified folder."))

include(joinpath(data_path, "config.jl"))

isfile(joinpath(data_path, read1_file)) || throw(AssertionError("Cannot find a valid file with the filename $read1_file. Please edit config.jl."))
isfile(joinpath(data_path, read2_file)) || throw(AssertionError("Cannot find a valid file with the filename $read2_file. Please edit config.jl."))
isfile(joinpath(data_path, annotation_file)) || throw(AssertionError("Cannot find a valid file with the filename $annotation_file. Please edit config.jl."))
isfile(joinpath(data_path, genome_file)) || throw(AssertionError("Cannot find a valid file with the filename $genome_file. Please edit config.jl."))

Pkg.activate(joinpath(@__DIR__, "ChimericAnalysis"))
Pkg.instantiate()

using ChimericAnalysis, BioSequences

genome = Genome(joinpath(data_path, genome_file))

features = Features(joinpath(data_path, annotation_file), [srna_type, cds_type, rrna_type, trna_type])
if autocomplete_utrs
    add5utrs!(features; utr_type=fiveutr_type, utr_length=autocomplete_utr_length)
    add3utrs!(features; utr_type=threeutr_type, utr_length=autocomplete_utr_length)
else
    merge!(features, Features(joinpath(data_path, annotation_file), [threeutr_type, fiveutr_type]))
end

if autocomplete_igrs
    addigrs!(features; igr_type=igr_type)
else
    merge!(features, Features(joinpath(data_path, annotation_file), [igr_type]))
end

show(features)

samplename_barcode = Dict(t[1] => LongDNA{4}(t[3]) for t in samplename_condition_antibarcode)
processing_path = mkpath(joinpath(data_path, "data_processing"))
fastqs = split_libs(joinpath(data_path, read1_file), joinpath(data_path, read2_file), samplename_barcode, processing_path)
bams = align_mem(fastqs, genome; bwa_bin=bwa_mem2_bin, sam_bin=samtools_bin)

results_path = mkpath(joinpath(data_path, "results"))
conditions = Dict(c => [i for (i, info) in enumerate(samplename_condition_antibarcode) if info[2]===c] for c in unique([t[2] for t in samplename_condition_antibarcode]))
chimeric_analysis(features, bams, results_path, conditions; filter_types=filter_types, min_distance=min_distance, prioritize_type=srna_type,
overwrite_type=igr_type, cds_type=cds_type, is_reverse_complement=is_reverse_complement, min_reads=min_reads, max_fdr=max_fdr,
include_read_identity=include_read_identity, include_singles=include_singles)
