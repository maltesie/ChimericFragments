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
        any(type(f) == cds_type for f in fs) || continue
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

function chimeric_analysis(features::Features, bams::SingleTypeFiles, results_path::String, conditions::Dict{String, Vector{Int}}, genome::Genome;
                            filter_types=["rRNA", "tRNA"],
                            min_distance=1000,
                            prioritize_type="sRNA",
                            min_prioritize_overlap=0.8,
                            max_bp_fdr=0.05,
                            overwrite_type="IGR",
                            max_ligation_distance=3,
                            is_reverse_complement=true,
                            is_paired_end=true,
                            check_interaction_distances=(30,0),
                            include_secondary_alignments=true,
                            include_alternative_alignments=false,
                            min_reads=5,
                            max_fisher_fdr=0.05,
                            fisher_exact_tail="right",
                            include_read_identity=true,
                            include_singles=true,
                            allow_self_chimeras=true,
                            bp_parameters=(4,5,0,7,8,3),
                            n_genome_samples=500000,
                            shift_weight=1.0,
                            keep_ints_without_ligation=true,
                            filter_name_queries=[])

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

        @info "Using $(summarize(features))"
        isdir(joinpath(results_path, "tables")) || mkpath(joinpath(results_path, "tables"))
        isdir(joinpath(results_path, "jld")) || mkpath(joinpath(results_path, "jld"))
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

            @info "Found $(length(interactions.multichimeras)) unique multichimeric arrangements."

            @info "Running statistical tests..."
            addpositions!(interactions, features)
            addpvalues!(interactions, genome, genome_model_ecdf; include_singles=include_singles, include_read_identity=include_read_identity,
                    fisher_exact_tail=fisher_exact_tail, check_interaction_distances=check_interaction_distances, bp_parameters=bp_parameters,
                    shift_weight=shift_weight)

            total_reads = sum(interactions.edges[!, :nb_ints])
            total_ints = nrow(interactions.edges)

            filterset = Set(findall(foldl(.|, occursin.(q, interactions.nodes.name) for q in filter_name_queries)))
            filter!([:src, :dst] => (x, y) -> !((x in filterset) || (y in filterset)), interactions.edges)

            above_min_reads = sum(interactions.edges[interactions.edges.nb_ints .>= min_reads, :nb_ints])
            above_min_ints = sum(interactions.edges.nb_ints .>= min_reads)
            filter!(:nb_ints => x -> x >= min_reads, interactions.edges)

            total_sig_reads = sum(interactions.edges[interactions.edges.fisher_fdr .<= max_fisher_fdr, :nb_ints])
            total_sig_ints = sum(interactions.edges.fisher_fdr .<= max_fisher_fdr)
            filter!(:fisher_fdr => x -> x <= max_fisher_fdr, interactions.edges)

            above_min_bp_reads = sum(interactions.edges[interactions.edges.bp_fdr .<= max_bp_fdr, :nb_ints])
            above_min_bp_ints = sum(interactions.edges.bp_fdr .<= max_bp_fdr)
            filter!(:bp_fdr => x -> (x <= max_bp_fdr) | (isnan(x) & keep_ints_without_ligation), interactions.edges)

            infotable = DataFrame(""=>["total interactions:", "annotation pairs:"], "total"=>[total_reads, total_ints], "reads>=$min_reads"=>[above_min_reads, above_min_ints],
                "& fdr<=$max_fisher_fdr"=>[total_sig_reads, total_sig_ints], "& bp_fdr<=$max_bp_fdr"=>[above_min_bp_reads, above_min_bp_ints])

            @info "interaction stats for condition $condition:\n" * DataFrames.pretty_table(String, infotable, nosubheader=true)

            correlation_matrix = corspearman(Matrix(interactions.edges[:, interactions.replicate_ids]))
            correlation_df = DataFrame(replicate_ids=interactions.replicate_ids)
            for (i,repid) in enumerate(interactions.replicate_ids)
                correlation_df[:, repid] = correlation_matrix[:, i]
            end
            @info "Correlation between interaction counts per replicate:\n" * DataFrames.pretty_table(String, correlation_df, nosubheader=true)

            @info "Saving..."

            odf = asdataframe(interactions; output=:edges, min_reads=min_reads, max_fisher_fdr=max_fisher_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "interactions_$(condition).csv"), odf)

            odf = asdataframe(interactions; output=:nodes, min_reads=min_reads, max_fisher_fdr=max_fisher_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "genes_$(condition).csv"), odf)

            odf = asdataframe(interactions; output=:multi)
            CSV.write(joinpath(results_path, "tables", "multi_$(condition).csv"), odf)

            write(joinpath(results_path, "jld", "$(condition).jld2"), interactions)

        end

        @info "Done."
    end
end