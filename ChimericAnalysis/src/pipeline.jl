function mergetypes(features::Features, cds_type::String, five_type::String, three_type::String, mergetype::String)
    # Compute merged annotations for anntoations of types cds_type, five_type and three_type with the same name
    types = (cds_type, five_type, three_type)
    merged_features = Interval{Annotation}[]

    # Collect all features with the types cds_type, five_type and three_type grouped by the annotation name
    feature_collector = Dict{String,Vector{Interval{Annotation}}}(name(feature)=>Interval{Annotation}[] for feature in features if type(feature) in types)
    for feature in features
        if type(feature) in types
            push!(feature_collector[name(feature)], feature)
        else
            push!(merged_features, feature)
        end
    end

    # Loop through grouped annotations and combine, if name, strand and reference match
    for (n, fs) in feature_collector
        s = strand(fs[1])
        ref = refname(fs[1])
        any((refname(f)!=ref) || (strand(f)!=s) for f in fs) &&
            throw(AssertionError("Features with the same name of types $types have to be on the same reference sequence and strand."))

        # Compute the new borders of the merged annotation
        left = minimum(leftposition(f) for f in fs)
        right = maximum(rightposition(f) for f in fs)

        # Can only merge if a cds_type is in the grouped annotations
        any(type(f) == cds_type for f in fs) || continue
        cds = first((s == STRAND_POS ? leftposition : rightposition)(f) for f in fs if type(f) == cds_type)

        # Save the position of the first nucleotide of the CDS for later use in visualisation
        ann = Dict("Name"=>n, "cds"=>isnothing(cds) ? "NA" : "$cds")
        push!(merged_features, Interval(ref, left, right, s, Annotation(mergetype, n, ann)))
    end
    return Features(merged_features)
end

function Base.write(fname::String, files::SingleTypeFiles)
    # Convert a set of CSV files into an XLSX file with all CSVs as sheets within.
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
                            min_self_chimera_distance=50,
                            bp_parameters=(4,5,0,7,8,3),
                            n_genome_samples=500000,
                            shift_weight=1.0,
                            keep_ints_without_ligation=true,
                            filter_name_queries=[],
                            join_pvalues_method=:stouffer,
                            min_mapping_quality=0)

    # Set up logger. It saves all @info and @warn messages and also outputs them to the terminal.
    filelogger = FormatLogger(joinpath(results_path, "analysis.log"); append=true) do io, args
        println(io, "[", args.level, "] ", args.message)
    end
    with_logger(TeeLogger(filelogger, ConsoleLogger())) do
        @info "Starting new analysis..."

        # Set up the AffineGapScoreModel used to compute the basepairing predictions. Parameters are defined in the config.jl
        scores = Dict((DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
                    (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
                    (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3])
        model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
            gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

        # The length of the sequences used for computing complementarity in the random model and the experiment. Parameters are defined in the config.jl
        seq_length = check_interaction_distances[1]-check_interaction_distances[2]+1

        # Compute a empirical density of the complementarity scores for randomly sampled pairs of sequences of fixed length.
        genome_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
            i1 % 2 == 0 ? genome.seq[i1:i1+seq_length-1] : reverse_complement(genome.seq[i1:i1+seq_length-1]),
            i2 % 2 == 0 ? reverse(genome.seq[i2:i2+seq_length-1]) : complement(genome.seq[i2:i2+seq_length-1]),
            model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genome.seq)-seq_length), (n_genome_samples,2)))]
        )

        # Set up folder structure for the results
        @info "Using $(summarize(features))"
        isdir(joinpath(results_path, "tables")) || mkpath(joinpath(results_path, "tables"))
        isdir(joinpath(results_path, "jld")) || mkpath(joinpath(results_path, "jld"))

        # Cycle through all conditions. The are defined in the config.jl as samplename_condition
        for (condition, r) in conditions

            replicate_ids = Vector{Symbol}()
            interactions = Interactions()
            GC.gc()

            @info "Collecting $(length(r)) samples for condition $condition:"

            # Reuse existing analysis results if present
            if isfile(joinpath(results_path, "jld", "$(condition).jld2"))
                @info "Found results files. Using existing JLD2 file for $condition..."
                interactions = Interactions(joinpath(results_path, "jld", "$(condition).jld2"))
            else
                # Cycle through all replicates in the current condition. The are defined in the config.jl as samplename_condition
                for (i, bam) in enumerate(bams[r])
                    replicate_id = Symbol("$(condition)_$i")
                    push!(replicate_ids, replicate_id)
                    @info "Replicate $replicate_id:"

                    @info "Reading $bam"
                    # Use the AlignedReads constructor from RNASeqTools to read alignments and sort them according to their position on the read.
                    alignments = AlignedReads(bam; include_secondary_alignments=include_secondary_alignments,
                                            include_alternative_alignments=include_alternative_alignments,
                                            is_reverse_complement=is_reverse_complement, min_mapping_quality=min_mapping_quality)

                    @info "Annotating alignments..."
                    # Use annotate! from RNASeqTools to efficiently and uniquely match the annotation with largest overlap to each alignment.
                    annotate!(alignments, features; prioritize_type=prioritize_type, min_prioritize_overlap=min_prioritize_overlap,
                                                    overwrite_type=overwrite_type)

                    @info "Building graph of interactions..."
                    # Call append! from interactions.jl to collect all interactions with their ligation points from the alignments.
                    append!(interactions, alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                        filter_types=filter_types, allow_self_chimeras=allow_self_chimeras, min_self_chimera_distance=min_self_chimera_distance,
                        is_paired_end=is_paired_end)

                    # Memory management
                    empty!(alignments)
                    GC.gc()
                end
            end

            length(interactions) == 0 && (@warn "No interactions found!"; continue)

            @info "Found $(length(interactions.multichimeras)) unique multichimeric arrangements."

            @info "Running statistical tests..."
            # Add infromation on relative positions of interactions within their annotations. Defined in interactions.jl
            addpositions!(interactions, features)

            # Filter data according to query strings which are looked up in the annotation names
            if !isempty(filter_name_queries)
                filterset = Set(findall(foldl(.|, occursin.(q, interactions.nodes.name) for q in filter_name_queries)))
                filter!([:src, :dst] => (x, y) -> !((x in filterset) || (y in filterset)), interactions.edges)
            end

            # Compute stats on significant interactions and read counts and filter data
            total_reads = sum(interactions.edges[!, :nb_ints])
            total_ints = nrow(interactions.edges)

            # Filter according to read count cutoff
            above_min_reads = sum(interactions.edges[interactions.edges.nb_ints .>= min_reads, :nb_ints])
            above_min_ints = sum(interactions.edges.nb_ints .>= min_reads)
            filter!(:nb_ints => x -> x >= min_reads, interactions.edges)

            # Do statistical tests. Defined in interactions.jl
            addpvalues!(interactions, genome, genome_model_ecdf; include_singles=include_singles, include_read_identity=include_read_identity,
                    fisher_exact_tail=fisher_exact_tail, check_interaction_distances=check_interaction_distances, bp_parameters=bp_parameters,
                    shift_weight=shift_weight, join_pvalues_method=join_pvalues_method)

            # Filter according to fishers exact test FDR cutoff
            total_sig_reads = sum(interactions.edges[interactions.edges.fisher_fdr .<= max_fisher_fdr, :nb_ints])
            total_sig_ints = sum(interactions.edges.fisher_fdr .<= max_fisher_fdr)
            filter!(:fisher_fdr => x -> x <= max_fisher_fdr, interactions.edges)

            # Filter according to complementarity FDR cutoff
            above_min_bp_reads = sum(interactions.edges[interactions.edges.bp_fdr .<= max_bp_fdr, :nb_ints])
            above_min_bp_ints = sum(interactions.edges.bp_fdr .<= max_bp_fdr)
            filter!(:bp_fdr => x -> (x <= max_bp_fdr) | (isnan(x) & keep_ints_without_ligation), interactions.edges)

            # Print pretty table
            #infotable = DataFrame(""=>["total interactions:", "annotation pairs:"], "total"=>[total_reads, total_ints], "reads>=$min_reads"=>[above_min_reads, above_min_ints],
            #    "& fdr<=$max_fisher_fdr"=>[total_sig_reads, total_sig_ints], "& bp_fdr<=$max_bp_fdr"=>[above_min_bp_reads, above_min_bp_ints])
            infotable = [
                "total interactions:" total_reads above_min_reads total_sig_reads above_min_bp_reads;
                "annotation pairs:" total_ints above_min_ints total_sig_ints above_min_bp_ints;
            ]

            @info "interaction stats for condition $condition:\n" * pretty_table(String, infotable,
                header=["", "total", "reads>=$min_reads", "& fdr<=$max_fisher_fdr", "& bp_fdr<=$max_bp_fdr"])

            # Check if interactions are left after filtering, compute correlation between interaction counts in replicates
            if nrow(interactions.edges) > 0
                correlation_matrix = hcat(String.(interactions.replicate_ids), cor(Matrix(interactions.edges[:, interactions.replicate_ids])))
                @info "Correlation between interaction counts per replicate:\n" * pretty_table(String, correlation_matrix,
                    header=vcat(["replicate_ids"], String.(interactions.replicate_ids)))
            else
                @info "Could not compute correlation between interaction counts. No interactions left after filtering!"
            end

            @info "Saving..."

            # Save interaction data in table form
            odf = asdataframe(interactions; output=:edges, min_reads=min_reads, max_fisher_fdr=max_fisher_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "interactions_$(condition).csv"), odf)

            # Save summary data per annotation in table form
            odf = asdataframe(interactions; output=:nodes, min_reads=min_reads, max_fisher_fdr=max_fisher_fdr, max_bp_fdr=max_bp_fdr)
            CSV.write(joinpath(results_path, "tables", "genes_$(condition).csv"), odf)

            # Save all found multi-chimeras into a table
            odf = asdataframe(interactions; output=:multi)
            CSV.write(joinpath(results_path, "tables", "multi_$(condition).csv"), odf)

            # Save the Interactions struct to JLD2 file
            write(joinpath(results_path, "jld", "$(condition).jld2"), interactions)

        end

        @info "Done."
    end
end