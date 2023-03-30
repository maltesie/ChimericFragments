function mergetypes(features::Features, cds_type::String, five_type::String, three_type::String, mergetype::String)
    types = (cds_type, five_type, three_type)
    merged_features = Interval{Annotation}[]
    feature_collector = Dict{String,Vector{Interval{Annotation}}}(name(feature)=>Interval{Annotation}[] for feature in features if type(feature) in types)
    for feature in features
        if type(feature) in types
            push!(feature_collector[name(feature)], feature)
        else
            push!(merged_features, feature)
        end
    end
    for (n, fs) in feature_collector
        s = strand(fs[1])
        ref = refname(fs[1])
        any((refname(f)!=ref) || (strand(f)!=s) for f in fs) &&
            throw(AssertionError("Features with the same name of types $types have to be on the same reference sequence and strand."))
        left = minimum(leftposition(f) for f in fs)
        right = maximum(rightposition(f) for f in fs)
        cds = first((s == STRAND_POS ? leftposition : rightposition)(f) for f in fs if type(f) == cds_type)
        ann = Dict("Name"=>n, "cds"=>isnothing(cds) ? "NA" : "$cds")
        push!(merged_features, Interval(ref, left, right, s, Annotation(mergetype, n, ann)))
    end
    return Features(merged_features)
end

function Base.write(fname::String, files::SingleTypeFiles)
    files.type in (".csv",) || throw(AssertionError("File type has to be .csv"))
    if files.type == ".csv"
        tables = Vector{Tuple{String,Vector{Any},Vector{String}}}()
        for file in files
            sheetname = basename(file)[1:end-length(files.type)]
            dataframe = DataFrame(CSV.File(file; stringtype=String))
            push!(tables,(sheetname,collect(eachcol(dataframe)), names(dataframe)))
        end
        XLSX.writetable(fname, tables; overwrite=true)
    end
end

score_bp(paln::PairwiseAlignmentResult, shift_weight::Float64) = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))

function chimeric_analysis(features::Features, bams::SingleTypeFiles, results_path::String, conditions::Dict{String, Vector{Int}}, genome::Genome;
                            filter_types=["rRNA", "tRNA"], min_distance=1000, prioritize_type="sRNA", min_prioritize_overlap=0.8, max_bp_fdr=0.05,
                            overwrite_type="IGR", max_ligation_distance=5, is_reverse_complement=true, is_paired_end=true, check_interaction_distances=(45,-10),
                            include_secondary_alignments=true, include_alternative_alignments=false, min_reads=5, max_fdr=0.05, fisher_exact_tail="right",
                            include_read_identity=true, include_singles=true, allow_self_chimeras=true, position_distribution_bins=50,
                            bp_parameters=(4,5,1,5,6,4), n_genome_samples=200000, plot_fdr_levels=[0.1,0.2,0.5], shift_weight=0.1)

    filelogger = FormatLogger(joinpath(results_path, "analysis.log"); append=true) do io, args
        println(io, "[", args.level, "] ", args.message)
    end
    with_logger(TeeLogger(filelogger, ConsoleLogger())) do
        @info "Starting new analysis..."

        scores = Dict((DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
                    (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
                    (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3])
        model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
            gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

        seq_length = check_interaction_distances[1]-check_interaction_distances[2]+1

        genome_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
            i1 % 2 == 0 ? genome.seq[i1:i1+seq_length-1] : reverse_complement(genome.seq[i1:i1+seq_length-1]),
            i2 % 2 == 0 ? reverse(genome.seq[i2:i2+seq_length-1]) : complement(genome.seq[i2:i2+seq_length-1]),
            model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genome.seq)-seq_length), (n_genome_samples,2)))]
        )

        randseq = randdnaseq(length(genome.seq))
        randseq_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
            view(randseq, i1:i1+seq_length-1),
            view(randseq, i2:i2+seq_length-1),
            model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genome.seq)-seq_length), (n_genome_samples,2)))]
        )

        @info "Using $(summarize(features))"
        isdir(joinpath(results_path, "tables")) || mkpath(joinpath(results_path, "tables"))
        isdir(joinpath(results_path, "jld")) || mkpath(joinpath(results_path, "jld"))
        isdir(joinpath(results_path, "plots")) || mkpath(joinpath(results_path, "plots"))
        for (condition, r) in conditions

            replicate_ids = Vector{Symbol}()
            interactions = Interactions()
            GC.gc()

            @info "Collecting $(length(r)) samples for condition $condition:"
            if isfile(joinpath(results_path, "jld", "$(condition).jld2"))
                @info "Found results files. Using existing JLD2 file for $condition..."
                interactions = Interactions(joinpath(results_path, "jld", "$(condition).jld2"))
            else
                for (i, bam) in enumerate(bams[r])
                    replicate_id = Symbol("$(condition)_$i")
                    push!(replicate_ids, replicate_id)
                    @info "Replicate $replicate_id:"
                    @info "Reading $bam"
                    alignments = AlignedReads(bam; include_secondary_alignments=include_secondary_alignments,
                                            include_alternative_alignments=include_alternative_alignments,
                                            is_reverse_complement=is_reverse_complement)
                    @info "Annotating alignments..."
                    annotate!(alignments, features; prioritize_type=prioritize_type, min_prioritize_overlap=min_prioritize_overlap,
                                                    overwrite_type=overwrite_type)
                    @info "Building graph of interactions..."
                    append!(interactions, alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                        filter_types=filter_types, allow_self_chimeras=allow_self_chimeras, is_paired_end=is_paired_end)
                    empty!(alignments)
                    GC.gc()
                end
            end

            length(interactions) == 0 && (@warn "No interactions found!"; continue)
            correlation_matrix = cor(Matrix(interactions.edges[:, interactions.replicate_ids]))
            correlation_df = DataFrame(replicate_ids=interactions.replicate_ids)
            for (i,repid) in enumerate(interactions.replicate_ids)
                correlation_df[:, repid] = correlation_matrix[:, i]
            end
            @info "Correlation between interaction counts:\n" * DataFrames.pretty_table(String, correlation_df, nosubheader=true)
            @info "Running statistical tests..."
            addpositions!(interactions, features)
            #ss = Vector{Set{Tuple{Int,Int}}}()
            #randpos = rand(1:(length(genome.seq)-50), (n_genome_samples,2))
            #for shift_weight in (0.0, 1.0, 1.5, 2.0)
            #    check_interaction_dists = check_interaction_distances
            #    bp_params = bp_parameters
            #    println("new test:")
            #    println(shift_weight)
            #    seq_length = check_interaction_dists[1]-check_interaction_dists[2]
            #    scores = Dict((DNA_A, DNA_T)=>bp_params[1], (DNA_T, DNA_A)=>bp_params[1],
            #        (DNA_C, DNA_G)=>bp_params[2], (DNA_G, DNA_C)=>bp_params[2],
            #        (DNA_G, DNA_T)=>bp_params[3], (DNA_T, DNA_G)=>bp_params[3])
            #    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_params[4], default_mismatch=-1*bp_params[4]);
            #        gap_open=-1*bp_params[5], gap_extend=-1*bp_params[6])
            #    genome_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
            #        i1 % 2 == 0 ? genome.seq[i1:i1+seq_length-1] : reverse_complement(genome.seq[i1:i1+seq_length-1]),
            #        i2 % 2 == 0 ? reverse(genome.seq[i2:i2+seq_length-1]) : complement(genome.seq[i2:i2+seq_length-1]),
            #        model), shift_weight) for (i1, i2) in eachrow(randpos)]
            #    )
            #    addpvalues!(interactions, genome, genome_model_ecdf; include_singles=include_singles, include_read_identity=include_read_identity,
            #        fisher_exact_tail=fisher_exact_tail, check_interaction_distances=check_interaction_dists, bp_parameters=bp_params, shift_weight=shift_weight)
            #    println(sum(interactions.edges.pred_fdr .<= 0.1))
            #    println(sum(interactions.edges.pred_fdr .<= 0.3))
            #    push!(ss, Set(collect(zip(interactions.edges.src[interactions.edges.pred_fdr .<= 0.1], interactions.edges.dst[interactions.edges.pred_fdr .<= 0.1]))))
            #    println()
            #end
            #println("overlaps between sets:")
            #println(length.(ss))
            #for r in eachrow([length(intersect(s1, s2)) for s1 in ss, s2 in ss])
            #    println(r)
            #end
            #println("normalized to first:")
            #for r in eachrow([round.(length(intersect(s1, s2))/length(s1); digits=2) for s1 in ss, s2 in ss])
            #    println(r)
            #end
            #println("normalized to min:")
            #for r in eachrow([round.(length(intersect(s1, s2))/min(length(s1), length(s2)); digits=2) for s1 in ss, s2 in ss])
            #    println(r)
            #end
            #supers = Set{Tuple{Int,Int}}()
            #for s in ss
            #    union!(supers, s)
            #end
            #println(length(supers))

            addpvalues!(interactions, genome, genome_model_ecdf; include_singles=include_singles, include_read_identity=include_read_identity,
                    fisher_exact_tail=fisher_exact_tail, check_interaction_distances=check_interaction_distances, bp_parameters=bp_parameters)

            total_reads = sum(interactions.edges[!, :nb_ints])
            total_ints = nrow(interactions.edges)

            above_min_reads = sum(interactions.edges[interactions.edges.nb_ints .>= min_reads, :nb_ints])
            above_min_ints = sum(interactions.edges.nb_ints .>= min_reads)

            total_sig_reads = sum(interactions.edges[interactions.edges.fdr .<= max_fdr, :nb_ints])
            total_sig_ints = sum(interactions.edges.fdr .<= max_fdr)

            both_reads = sum(interactions.edges[(interactions.edges.fdr .<= max_fdr) .& (interactions.edges.nb_ints .>= min_reads), :nb_ints])
            both_ints = sum((interactions.edges.fdr .<= max_fdr) .& (interactions.edges.nb_ints .>= min_reads))

            above_min_bp_reads = sum(interactions.edges[interactions.edges.pred_fdr .<= max_bp_fdr, :nb_ints])
            above_min_bp_ints = sum(interactions.edges.pred_fdr .<= max_bp_fdr)

            infotable = DataFrame(""=>["total interactions:", "annotation pairs:"], "total"=>[total_reads, total_ints], "reads>=$min_reads"=>[above_min_reads, above_min_ints],
                "fdr<=$max_fdr"=>[total_sig_reads, total_sig_ints], "both"=>[both_reads, both_ints], "bp_fdr<=$max_bp_fdr"=>[above_min_bp_reads, above_min_bp_ints])

            @info "interaction stats for condition $condition:\n" * DataFrames.pretty_table(String, infotable, nosubheader=true)
            @info "Saving tables and plots..."
            odf = asdataframe(interactions; output=:edges, min_reads=min_reads, max_fdr=max_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "interactions_$(condition).csv"), odf)
            odf = asdataframe(interactions; output=:stats, min_reads=min_reads, max_fdr=max_fdr, max_bp_fdr=max_bp_fdr, hist_bins=position_distribution_bins)
            CSV.write(joinpath(results_path, "tables", "ligation_points_$(condition).csv"), odf)
            odf = asdataframe(interactions; output=:nodes, min_reads=min_reads, max_fdr=max_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "genes_$(condition).csv"), odf)

            write(joinpath(results_path, "jld", "$(condition).jld2"), interactions)

            p = bp_score_dist_plot(interactions, genome_model_ecdf, randseq_model_ecdf, plot_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_bp_scores_dist.png"), p)

            (p2_1, p2_2, p2_3), (p3_1, p3_2, p3_3) = bp_clipping_dist_plots(interactions, check_interaction_distances, plot_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_clippings_bp.png"), p2_1)
            save(joinpath(results_path, "plots", "$(condition)_clippings_fisher.png"), p3_1)
            save(joinpath(results_path, "plots", "$(condition)_clippings_bp_ratios.png"), p2_2)
            save(joinpath(results_path, "plots", "$(condition)_clippings_fisher_ratios.png"), p3_2)
            save(joinpath(results_path, "plots", "$(condition)_clippings_bp_diff.png"), p2_3)
            save(joinpath(results_path, "plots", "$(condition)_clippings_fisher_diff.png"), p3_3)

            p4, p5 = interaction_distribution_plots(interactions, plot_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_count_dist.png"), p4)
            save(joinpath(results_path, "plots", "$(condition)_odds_ratio_dist.png"), p5)

            p6 = node_distribution_plot(interactions, plot_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_degree_dist.png"), p6)

            p7 = annotation_type_heatmap(interactions, plot_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_annotation_type_heatmap.png"), p7)
        end

        @info "Done."
    end
end