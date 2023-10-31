# Data structure to hold all information on all interactions including ligation points and basepairing predictions.
struct Interactions
    # Table connecting an index with an annotation and summary information on interactions
    nodes::DataFrame
    # Table containing all interactions with interaction counts split by replicate
    edges::DataFrame
    # Dictionary containing all interacting pairs as keys connected with information on ligation points and interaction sites
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    # Dictionary containing all ligation points as key connected with additional information such as p-value, position of the complementarity region
    # within the sequences and complementarity score.
    bpstats::Dict{Tuple{Int, Int, Int, Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int64}}
    # Dictionary counting all multi-chimeric arrangements
    multichimeras::Dict{Vector{Int}, Int}
    # Replicate IDs from the config get stored here for access of each replicates' counts in the edges table
    replicate_ids::Vector{Symbol}
    # Counts for the classification of reads into single, chimeric, self-chimeric, multi-chimeric and unclassified per replicate.
    counts::Dict{Symbol,Vector{Int}}
end

# Empty constructor
function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_selfchimeric=>Int[], :nb_unclassified=>Int[],
        :nb_ints=>Int[], :nb_ints_src=>Int[], :nb_ints_dst=>Int[], :nb_partners=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :meanlen1=>Float64[], :meanlen2=>Float64[], :nms1=>Float64[], :nms2=>Float64[])
    edgestats = Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}()
    bpstats = Dict{Tuple{Int, Int, Int, Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int64}}()
    multichimeras = Dict{Vector{Int}, Int}()
    counts = Dict{Symbol,Vector{Int}}()
    Interactions(nodes, edges, edgestats, bpstats, multichimeras, Symbol[], counts)
end

# Constructor for the first AlignedReads to be added to an empty Interactions struct
Interactions(alignments::AlignedReads; replicate_id=:first, min_distance=1000, max_ligation_distance=5, filter_types=[], allow_self_chimeras=false) =
    append!(Interactions(), alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                filter_types=filter_types, allow_self_chimeras=allow_self_chimeras)

# Count number of interaction in interactions
Base.length(interactions::Interactions) = nrow(interactions.edges)

# Deep memory clean. Julias garbage collector sometimes needs this to no reserve too much memory.
function Base.empty!(interactions::Interactions)
    empty!(interactions.nodes)
    empty!(interactions.edges)
    empty!(interactions.edgestats)
    empty!(interactions.bpstats)
    empty!(interactions.multichimeras)
    empty!(interactions.replicate_ids)
    empty!(interactions.counts)
end


# Method of write function which saves the Interactions struct in a jld2 file.
function Base.write(filepath::String, interactions::Interactions)
    if !endswith(filepath, ".jld2")
        throw(ArgumentError("Append '.jld2' to filepath"))
    else
        save(filepath, "interactions", interactions)
    end
end

# Open a jld2 file and parse it into an Interactions struct
Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

# Main function for adding a new replicate into an Interactions struct.
function Base.append!(interactions::Interactions, alignments::AlignedReads, replicate_id::Symbol;
    min_distance=1000, max_ligation_distance=5, filter_types=["rRNA", "tRNA"], allow_self_chimeras=false, is_paired_end=true)

    # If the supplied replicate_id is not known, instantiate new counter and add replicate_id to list of known ids
    if !(String(replicate_id) in interactions.replicate_ids)
        interactions.edges[:, replicate_id] = zeros(Int, nrow(interactions.edges))
        push!(interactions.replicate_ids, replicate_id)
    end

    counts = zeros(Int, 6) #single, unclassified, selfchimeric, exclude, chimeric, multi

    # Translation dictionary to link hashs of annotation names with corresponding index in the interactions table
    trans = Dict{UInt, Int}(interactions.nodes.hash[i]=>i for i in 1:nrow(interactions.nodes))
    edgestats = interactions.edgestats

    # Create Matrix from Dataframe for efficient manipulation
    node_ints = Matrix{Int}(interactions.nodes[:, [:nb_single, :nb_unclassified, :nb_selfchimeric, :nb_ints, :nb_ints_src, :nb_ints_dst]])

    # Preallocate some memory for new annotations and reuse info from already processed replicates.
    node_ids = zeros(Int, length(interactions.nodes.hash))
    node_hashs = Vector{UInt}(interactions.nodes.hash)
    nb_nodes_start = nrow(interactions.nodes)

    # Preallocate memory for counting interactions in the current replicate and copy counts from already processed replicates
    edge_ints = hcat(Matrix{Int}(interactions.edges[:, [:src, :dst, :nb_ints, :nb_multi]]), zeros(Int, nrow(interactions.edges)))
    edge_floats = Matrix{Float64}(interactions.edges[:, [:meanlen1, :meanlen2, :nms1, :nms2]])

    # Preallocate a MergedAlignedRead (defined in mergedreads.jl) to reuse for each read
    mergedread = MergedAlignedRead(alignments)

    # Define hash variable
    h = UInt(1)

    # Loop through all reads, the alignment variable is a AlignedRead from RNASeqTools containing all alignments from one read or read pair
    for alignment in alignments
        
         # Skip read if it contains no annotated fragments
        if !(hasannotation(alignment))
            counts[2] += 1
            continue
        end

        # Compute merged alignments, combining overlapping information from reads 1 and 2.
        mergeparts!(mergedread, alignment)

        # Skip read if it contains only alignments annotated with a type defined in filter_types
        if !isempty(filter_types) && all(alignments.antypes[first(pair)] in filter_types for pair in mergedread.pindexpairs)
            counts[4] += 1
            continue
        end

        # Classify read into chimeric and possibly multi-chimeric
        is_chimeric = ischimeric(mergedread; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras)
        counts[5] += is_chimeric
        is_multi = length(mergedread)>2

        # Loop through all merged alignments to add hashs to the translation dictionary for indexing
        for (i1,_) in mergedread.pindexpairs

            # Skip if type of annotation is in filter_types
            !isempty(filter_types) && (alignments.antypes[i1] in filter_types) && continue
            
            # Compute hash based on name and type of the annotation
            h = myhash(alignments, i1)
            if !(h in keys(trans))
                idx = length(trans) + 1

                # Connect hash with index in nodes table
                trans[h] = idx

                # allocate memory in chunks
                if idx > length(node_hashs)
                    resize!(node_hashs, length(node_hashs) + 100000)
                    resize!(node_ids, length(node_hashs) + 100000)
                    node_ints = vcat(node_ints, zeros(Int, 100000, 6))
                end

                node_ids[idx] = i1
                node_hashs[idx] = h
            end
        end

        # Classify read into self-chimeric xor single xor leave unclassified and count accordingly
        if isselfchimeric(mergedread; min_distance=min_distance)
            node_ints[trans[h], 3] += 1
            counts[3] += 1
        elseif !is_chimeric
            if issingle(mergedread; is_paired_end=is_paired_end)
                node_ints[trans[h], 1] += 1
                counts[1] += 1
            else
                node_ints[trans[h], 2] += 1
                counts[2] += 1
            end
        end

        # Add multi-chimeras to dictionary using the indices in the node table in the correct order on the read(s) as a key.
        if is_multi
            all_together = [trans[myhash(alignments, i1)] 
                for (i1,_) in mergedread.pindexpairs if myhash(alignments, i1) in keys(trans)]
            if length(all_together) > 2 
                counts[6] += 1
                if all_together in keys(interactions.multichimeras)
                    interactions.multichimeras[all_together] += 1
                else
                    interactions.multichimeras[all_together] = 1
                end
            end
        end

        # Loop through all possible combinations of aligned parts on reads while keeping the order on the read. pair1 is upstream of pair2 on the read.
        for (i, pair1) in enumerate(mergedread.pindexpairs), pair2 in (@view mergedread.pindexpairs[i+1:end])

            # Skip if pair is not chimeric. This is used to deal with self-chimeras.
            ischimeric(mergedread, pair1, pair2; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras) || continue

             # Skip if type of annotation is in filter_types
            (!isempty(filter_types) && (alignments.antypes[first(pair1)] in filter_types || alignments.antypes[first(pair2)] in filter_types)) && continue
            
            # Get indices from hashs of the annotations of both alignments
            a, b = trans[myhash(alignments, first(pair1))], trans[myhash(alignments, first(pair2))]

            # Count summary for each annotation with info on its order in the pair.
            node_ints[a, 4] += 1
            node_ints[a, 5] += 1
            node_ints[b, 4] += 1
            node_ints[b, 6] += 1

            # Allocate new ligation points counter if interaction was not seen before
            if !((a,b) in keys(edgestats))
                idx = length(edgestats)+1

                # Instantiate new ligation points counter
                edgestats[(a,b)] = (idx, Dict{Tuple{Int, Int},Int}(), Dict{Tuple{Int, Int},Int}())

                # Preallocate space for new interactions in chunks
                if idx > size(edge_ints)[1]
                    edge_ints = vcat(edge_ints, zeros(Int, 100000, 5))
                    edge_floats = vcat(edge_floats, zeros(Float64, 100000, 4))
                end

                # Connect interaction with corresponding indices in node table
                edge_ints[idx, 1] = a
                edge_ints[idx, 2] = b
            end

            # Fetch references to ligation point and interaction point counters
            (iindex, intcounter, ligationcounter) = edgestats[(a, b)]

            # Set global interaction count, multi-chimera count and current replicate count
            edge_ints[iindex, 3] += 1
            edge_ints[iindex, 4] += is_multi
            edge_ints[iindex, 5] += 1

            # Compute ligation point or interaction point. Its the same procedure for both, later its decided if this is saved as a ligation point or not.
            # Both are defined as a set of two coordinates on the genome. First, the corrdinate of the most upstream nucleotide of the downstream fragment
            # on the read. And second the most downstream coordinate of the upstream fragment.
            leftpos = alignments.strands[last(pair1)] === STRAND_NEG ? alignments.leftpos[last(pair1)] : alignments.rightpos[last(pair1)]
            rightpos = alignments.strands[first(pair2)] === STRAND_NEG ? alignments.rightpos[first(pair2)] : alignments.leftpos[first(pair2)]

            # Decide if the coordinates are saved as ligation point or interaction point. Its a ligation point if both fragments come from the same read and
            # have a maximum of max_ligation_distance nucleotides between them.
            counter = hasligationpoint(mergedread, pair1, pair2; max_distance=max_ligation_distance) ? ligationcounter : intcounter
            pospair = (leftpos, rightpos)
            pospair in keys(counter) ? (counter[pospair]+=1) : (counter[pospair]=1)

            # Save all averages by computing them as a running average using the known number of parts already in the average (unneccessarily complicated).
            # Saves average number of missmatches in the alignment and length of the alignment.
            for (i,v) in enumerate((reflen(alignments, pair1), reflen(alignments, pair2),
                            max(alignments.nms[first(pair1)], alignments.nms[last(pair1)]),
                            max(alignments.nms[first(pair2)], alignments.nms[last(pair2)])))
                edge_floats[iindex, i] += (v - edge_floats[iindex, i]) / edge_ints[iindex, 5]
            end
        end
    end

    # Save current replicates classification counts to the supplied Interactions struct.
    interactions.counts[replicate_id] = counts

    # Efficiently integrate new annotations from AlignedReads alignments container into the Interactions struct
    resize!(interactions.nodes, length(trans))
    nr = (nb_nodes_start+1):length(trans)
    view(interactions.nodes.name, nr) .= view(alignments.annames, view(node_ids, nr))
    view(interactions.nodes.type, nr) .= view(alignments.antypes, view(node_ids, nr))
    view(interactions.nodes.ref, nr) .= view(alignments.refnames, view(node_ids, nr))
    view(interactions.nodes.strand, nr) .= view(alignments.strands, view(node_ids, nr))
    view(interactions.nodes.hash, nr) .= view(node_hashs, nr)

    # Integrate summary of classification counts into nodes table
    nr = 1:length(trans)
    for (i, n) in enumerate((:nb_single, :nb_selfchimeric, :nb_unclassified, :nb_ints, :nb_ints_src, :nb_ints_dst))
        interactions.nodes[!, n] .= view(node_ints, nr, i)
    end

    # Efficiently integrate interaction counts into edges table
    er = 1:length(edgestats)
    nbefore = nrow(interactions.edges)
    resize!(interactions.edges, length(edgestats))
    for rid in interactions.replicate_ids
        view(interactions.edges[!, rid], nbefore+1:length(edgestats)) .= zeros(Int, length(edgestats)-nbefore)
    end
    for (i, n) in enumerate((:src, :dst, :nb_ints, :nb_multi))
        interactions.edges[!, n] .= view(edge_ints, er, i)
    end
    for (i, n) in enumerate((:meanlen1, :meanlen2, :nms1, :nms2))
        interactions.edges[!, n] .= view(edge_floats, er, i)
    end
    interactions.edges[!, replicate_id] .= view(edge_ints, er, 5)

    # Create info table for pretty output to terminal and log file
    infotable = DataFrame(
        :total=>[nread(alignments)], :chimeric=>[counts[5]], :multi=>[counts[6]], :self=>[counts[3]],
        :single=>[counts[1]], :filtered=>[counts[4]], :no_class=>[counts[2]],
    )
    @info "Classification of reads:\n" * DataFrames.pretty_table(String, infotable, nosubheader=true)
    return interactions
end

# Compute the complementarity score as a combination between the alignment score from the AffineGapScoreModel and the weighted shift between
# complementarity regions towards eachother in the local coordinates of the fixed-length sequences used to compute the complementarity.
score_bp(paln::PairwiseAlignmentResult, shift_weight::Float64) = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))

# Perform statistical tests (Fisher exact and custom complementarity test)
function addpvalues!(interactions::Interactions, genome::Genome, random_model_ecdf::ECDF; fisher_exact_tail="right", include_read_identity=true,
                        include_singles=true, check_interaction_distances=(30,0), bp_parameters=(4,5,0,7,8,3), shift_weight=0.5, join_pvalues_method=:stouffer)

    # Decide, how to populate contingency table for Fishers exact test. include_read_identity toggles, if the order in which the pair was found in on the read
    # is taken into account.
    if include_read_identity
        # Compute one p-value for every interaction. ints_between is the number of interactions for a pair of annotations a => b
        ints_between = interactions.edges[!, :nb_ints]
        # other_source is the number of interactions a is involved in as partner 1 (other than ints_between)
        other_source = interactions.nodes[interactions.edges[!, :src], :nb_ints_src] .- ints_between
        # other_target is the number of interactions b is involved in as partner 2 (other than ints_between)
        other_target = interactions.nodes[interactions.edges[!, :dst], :nb_ints_dst] .- ints_between
    else
        # Interactions a => b and b => a share a p-value, both orientations lead to the same results here.
        check_dict = Dict((s,d)=>n for (s,d,n) in eachrow(interactions.edges[!, [:src, :dst, :nb_ints]]))
        # ints_between is the number of interactions for a pair of annotations in each orientation (a => b and b => a)
        ints_between = [(d,s) in keys(check_dict) ? check_dict[(d,s)]+ints : ints for (s,d,ints) in eachrow(interactions.edges[!, [:src, :dst, :nb_ints]])]
        # other_source is the number of interactions a is involved in (other than ints_between)
        other_source = interactions.nodes[interactions.edges[!, :src], :nb_ints] .- ints_between
        # other_target is the number of interactions b is involved in (other than ints_between)
        other_target = interactions.nodes[interactions.edges[!, :dst], :nb_ints] .- ints_between
    end

    # Decide, how to populate contingency table for Fishers exact test. include_singles toggles, if single counts of the partners alone are taken into account.
    if include_singles
        other_source .+= interactions.nodes[interactions.edges[!, :src], :nb_single]
        other_source .+= interactions.nodes[interactions.edges[!, :src], :nb_selfchimeric]
        other_target .+= interactions.nodes[interactions.edges[!, :dst], :nb_single]
        other_target .+= interactions.nodes[interactions.edges[!, :dst], :nb_selfchimeric]
    end

    # Number of interactions not involving a or b
    total_other = sum(interactions.edges[!, :nb_ints]) .- ints_between .- other_source .- other_target .+
        (include_singles ? sum(interactions.nodes[!, :nb_single] .+ interactions.nodes[!, :nb_selfchimeric]) : 0)

    # Broadcast all tests and save odds ratio in the contingency table and the p-values and corresponding FDR
    tests = FisherExactTest.(ints_between, other_target, other_source, total_other)
    odds_ratio = [t.ω for t in tests]
    pvalues_fisher = pvalue.(tests; tail=Symbol(fisher_exact_tail))
    adjp_fisher = adjust(PValues(pvalues_fisher), BenjaminiHochberg())
    interactions.edges[:, :odds_ratio] = odds_ratio
    interactions.edges[:, :fisher_pvalue] = pvalues_fisher
    interactions.edges[:, :fisher_fdr] = adjp_fisher

    # preallocate memory for complementarity p-values and FDR
    interactions.edges[:, :bp_pvalue] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :bp_fdr] = fill(NaN, nrow(interactions.edges))

    # Setup AffineGapScoreModel for computing the complementarity. Parameters are defined in config.jl
    scores = Dict((DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
                    (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
                    (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3])
    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
        gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

    # Precompute complement and reversed genome for efficiency.
    complement_genome = Genome(complement(genome.seq), genome.chroms)
    reverse_genome = Genome(copy(genome.seq), genome.chroms)
    for (_, seq) in reverse_genome
        reverse!(seq)
    end
    reverse_complement_genome = Genome(complement(genome.seq), genome.chroms)
    for (_, seq) in reverse_complement_genome
        reverse!(seq)
    end

    # Loop through all interactions
    for (c, edge_row) in enumerate(eachrow(interactions.edges))

        # check if any ligation point or interaction point is in the data to prevent key errors
        if (edge_row.src, edge_row.dst) in keys(interactions.edgestats)

            # Loop through all ligation points (third entry of an edgestats record)
            for ((i1::Int, i2::Int), ligation_count::Int) in interactions.edgestats[(edge_row.src, edge_row.dst)][3]

                # Collect info about the interaction from the nodes table
                strand1, strand2 = interactions.nodes[edge_row[:src], :strand], interactions.nodes[edge_row[:dst], :strand]
                ref1, ref2 = interactions.nodes[edge_row[:src], :ref], interactions.nodes[edge_row[:dst], :ref]
                l1, l2 = length(genome.chroms[ref1]), length(genome.chroms[ref2])

                # Setup reference to fixed-length sequence of partner 1 from the corresponding genome (reversed or complemented)
                s1 = if strand1=='+'
                    view(genome[ref1], clamp(i1-check_interaction_distances[1]+1, 1, l1):clamp(i1-check_interaction_distances[2], 1, l1))
                else
                    c1, c2 = l1+1-clamp(i1+check_interaction_distances[2], 1, l1), l1+1-clamp(i1+check_interaction_distances[1]-1, 1, l1)
                    view(reverse_complement_genome[ref1], c2:c1)
                end

                # Setup reference to reversed fixed-length sequence of partner 2 from the corresponding genome (reversed or complemented)
                s2 = if strand2=='-'
                    view(complement_genome[ref2], clamp(i2-check_interaction_distances[1]+1, 1, l2):clamp(i2-check_interaction_distances[2], 1, l2))
                else
                    c1, c2 = l2+1-clamp(i2+check_interaction_distances[2], 1, l2), l2+1-clamp(i2+check_interaction_distances[1]-1, 1, l2)
                    view(reverse_genome[ref2], c2:c1)
                end

                # Compute local alignment with AffineGapScoreModel model to find best complementarity region, parameters get defined in config.jl
                paln = pairalign(LocalAlignment(), s1, s2, model)

                # Compute score and p-value from empirical cumulative density of random scores and save all into the Interactions struct.
                sco = score_bp(paln, shift_weight)
                thisp = 1-random_model_ecdf(sco)
                interactions.bpstats[(edge_row.src, i1, edge_row.dst, i2)] = (thisp, paln.aln.a.aln.anchors[1].seqpos + 1, paln.aln.a.aln.anchors[end].seqpos,
                                                        paln.aln.a.aln.anchors[1].refpos + 1, paln.aln.a.aln.anchors[end].refpos, sco, ligation_count)

            end

            # Combine p-values for each interaction from all corresponding ligation points. Possible methods are: Fisher, Stouffer and the minimum of FDR values
            if join_pvalues_method in (:fisher, :fdr)
                # Preallocate memory for p-value collection.
                pvs = zeros(sum(interactions.bpstats[(edge_row.src, p[1], edge_row.dst, p[2])][7] for p in keys(interactions.edgestats[(edge_row.src, edge_row.dst)][3])
                    if (edge_row.src, p[1], edge_row.dst, p[2]) in keys(interactions.bpstats); init=0))
                edge_row.bp_pvalue = if length(pvs) > 0
                    current_i = 1
                    # Populate preallocated p-values array by saving one p-value multiple times if the corresponding ligation point was sampled multiple times.
                    for p in keys(interactions.edgestats[(edge_row.src, edge_row.dst)][3])
                        bpkey = (edge_row.src, p[1], edge_row.dst, p[2])
                        if bpkey in keys(interactions.bpstats)
                            pvs[current_i:(current_i+interactions.bpstats[bpkey][7]-1)] .= interactions.bpstats[bpkey][1]
                        end
                        current_i += interactions.bpstats[bpkey][7]
                    end
                    # Compute combined p-value with Fishers method or as minimum of FDR values
                    join_pvalues_method == :fisher ? MultipleTesting.combine(PValues(pvs), Fisher()) : minimum(adjust(PValues(pvs), BenjaminiHochberg()))
                else
                    # Save a NaN if no ligation points were found for this interaction.
                    NaN
                end
            elseif join_pvalues_method == :stouffer
                # Stouffers method can weight the individual p-values by the number, the ligation point was sampled, so they can be collected directly.
                pvs = [interactions.bpstats[(edge_row.src, p[1], edge_row.dst, p[2])][1] for p in keys(interactions.edgestats[(edge_row.src, edge_row.dst)][3])
                    if (edge_row.src, p[1], edge_row.dst, p[2]) in keys(interactions.bpstats)]
                ws = Float64[interactions.bpstats[(edge_row.src, p[1], edge_row.dst, p[2])][7] for p in keys(interactions.edgestats[(edge_row.src, edge_row.dst)][3])
                    if (edge_row.src, p[1], edge_row.dst, p[2]) in keys(interactions.bpstats)]
                # Compute combined p-value by Stouffers method or return NaN if no ligation point was found.
                edge_row.bp_pvalue = length(pvs) > 0 ? MultipleTesting.combine(PValues(pvs), ws, Stouffer()) : NaN
            end
        end
    end

    # Compute FDR values again for the whole dataset.
    nan_index = .!isnan.(interactions.edges.bp_pvalue)
    any(nan_index) && (interactions.edges.bp_fdr[nan_index] = adjust(PValues(interactions.edges.bp_pvalue[nan_index]), BenjaminiHochberg()))
    return interactions
end

# Extract +1 coordinate from CDS from merged annotations (5UTR, 3UTR and CDS)
cdsposition(feature::Interval{Annotation}) = hasparam(feature, "cds") ? param(feature, "cds", Int) : 0

# Add coordinates of annotations to nodes table
function addpositions!(interactions::Interactions, features::Features)
    tus = Dict(hash(name(feature), hash(type(feature)))=>(leftposition(feature), rightposition(feature), cdsposition(feature)) for feature in features)
    interactions.nodes[:, :left] = Vector{Int}(undef, nrow(interactions.nodes))
    interactions.nodes[:, :right] = Vector{Int}(undef, nrow(interactions.nodes))
    interactions.nodes[:, :cds] = Vector{Int}(undef, nrow(interactions.nodes))
    for nodes_row in eachrow(interactions.nodes)
        nodes_row[[:left, :right, :cds]] = tus[nodes_row[:hash]]
    end
    return interactions
end

# Compute normalized histogram of ligation points or interaction points over the length of the anntoation.
function histo(ints::Dict{Tuple{Int,Int},Int}, mi::Int, ma::Int, nbins::Int, l::Bool)
    h = zeros(Int, nbins)
    dbin = (ma-mi+1)/nbins
    for i in mi:ma
        for (k,c) in ints
            if i == (l ? first(k) : last(k))
                h[Int(floor((i-mi)/dbin)) + 1] += c
            end
        end
    end
    return h
end
# Collect additional information to output interactions or summary information for all annotations.
function asdataframe(interactions::Interactions; output=:edges, min_reads=5, max_fisher_fdr=0.05, max_bp_fdr=0.05, hist_bins=100)
    # Filter output data by number of reads per interaction or Fisher FDR or complementarity FDR
    filter_index = (interactions.edges[!, :nb_ints] .>= min_reads) .& (interactions.edges[!, :fisher_fdr] .<= max_fisher_fdr) .&
                        ((interactions.edges[!, :bp_fdr] .<= max_bp_fdr) .| isnan.(interactions.edges.bp_fdr))
    out_df = interactions.edges[filter_index, :]

    # To output interactions table, add annotation details for each partner and format numbers for human readability.
    if output === :edges
        out_df[!, :meanlen1] = Int.(round.(out_df[!, :meanlen1]))
        out_df[!, :meanlen2] = Int.(round.(out_df[!, :meanlen2]))
        out_df[!, :nms1] = round.(out_df[!, :nms1], digits=4)
        out_df[!, :nms2] = round.(out_df[!, :nms2], digits=4)
        out_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        out_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        out_df[:, :ref1] = interactions.nodes[out_df[!,:src], :ref]
        out_df[:, :ref2] = interactions.nodes[out_df[!,:dst], :ref]
        out_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        out_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        out_df[:, :strand1] = interactions.nodes[out_df[!,:src], :strand]
        out_df[:, :strand2] = interactions.nodes[out_df[!,:dst], :strand]
        out_df[:, :left1] = interactions.nodes[out_df[!,:src], :left]
        out_df[:, :left2] = interactions.nodes[out_df[!,:dst], :left]
        out_df[:, :right1] = interactions.nodes[out_df[!,:src], :right]
        out_df[:, :right2] = interactions.nodes[out_df[!,:dst], :right]
        out_df[:, :in_libs] = sum(eachcol(out_df[!, interactions.replicate_ids] .!= 0))
        # Set column order
        out_columns = [:name1, :type1, :ref1, :strand1,:left1, :right1, :name2, :type2, :ref2, :strand2, :left2, :right2, :nb_ints, :nb_multi, :in_libs,
        :fisher_pvalue, :fisher_fdr, :odds_ratio, :bp_pvalue, :bp_fdr, :meanlen1, :nms1, :meanlen2, :nms2]
        return sort!(out_df[!, out_columns], :nb_ints; rev=true)

    # For output of summary information per annotation, compute sum of reads per interaction and number of partners
    elseif output === :nodes
        out_nodes = copy(interactions.nodes)
        for (i,row) in enumerate(eachrow(out_nodes))
            row[:nb_ints] = sum(out_df[(out_df.src .== i) .| (out_df.dst .== i), :nb_ints])
            partners = union!(Set(out_df[out_df.src .== i, :dst]), Set(out_df[out_df.dst .== i, :src]))
            row[:nb_partners] = length(partners)
        end
        # Set column order
        return sort!(out_nodes[!, [:name, :type, :ref, :nb_single, :nb_selfchimeric, :nb_unclassified, :nb_ints, :nb_partners]], :nb_single; rev=true)

    # For output of multi-chimera table, translate indices into annotations and create a table with number of columns equal to the maximum number of
    # partners found on a single read or read pair
    elseif output === :multi
        odf = DataFrame()
        if length(interactions.multichimeras) > 0
            # Reverse sort by number of times, this exact arrangement was found
            sorted_multis = sort(collect(interactions.multichimeras), by=x->x[2], rev=true)

            # Compute number of columns
            max_nb_partners = maximum(length(m[1]) for m in sorted_multis)

            # Make columns an enumerate them
            for i in 1:max_nb_partners
                odf[:, "name_$i"] = String[]
                odf[:, "type_$i"] = String[]
            end
            # Number of times, this arrangement was found
            odf[:, :count] = Int[]
            # Number of unique annotations involved in this chimera
            odf[:, :nb_partners] = Int[]
            current_multi = Vector{String}(undef, 2*max_nb_partners)
            for (multichimera, count) in sort(collect(interactions.multichimeras), by=x->x[2], rev=true)
                i = 0
                # Populate as many columns as the multi-chimera has partners (can be multiple instances of the same annotation)
                for mc in multichimera
                    current_multi[i+=1] = interactions.nodes.name[mc]
                    current_multi[i+=1] = interactions.nodes.type[mc]
                end
                # Add empty string for remaining columns
                current_multi[(i+1):(2*max_nb_partners)] .= ""
                # Push row (including chimera count and number of unique annotation involved in the chimera) to table
                push!(odf, (current_multi..., count, length(unique(multichimera))))
            end
        end
        return odf

    # To output a table with histogram of ligation points along normalzed length of all annotations
    elseif output === :stats
        edgestats = interactions.edgestats

        # Preallocate Matrix for histogram data for all annotations
        statsmatrix = zeros(nrow(out_df), 4*hist_bins)

        # Set up columns according to number of bins in histogram
        stats_df = DataFrame(statsmatrix, [
            ["$(i)_ints1" for i in 1:hist_bins]...,
            ["$(i)_ints2" for i in 1:hist_bins]...,
            ["$(i)_lig1" for i in 1:hist_bins]...,
            ["$(i)_lig2" for i in 1:hist_bins]...
            ])

        # Add annotation details for each partner
        stats_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        stats_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        stats_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        stats_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        stats_df[:, :strand1] = interactions.nodes[out_df[!,:src], :strand]
        stats_df[:, :strand2] = interactions.nodes[out_df[!,:dst], :strand]
        stats_df[:, :left1] = interactions.nodes[out_df[!,:src], :left]
        stats_df[:, :left2] = interactions.nodes[out_df[!,:dst], :left]
        stats_df[:, :right1] = interactions.nodes[out_df[!,:src], :right]
        stats_df[:, :right2] = interactions.nodes[out_df[!,:dst], :right]

        # Populate matrix with histogram data for every annotation.
        for (i, (a,b, mi1, ma1, mi2, ma2)) in enumerate(zip(out_df[!,:src], out_df[!,:dst], out_df[!, :left1], out_df[!, :right1], out_df[!, :left2], out_df[!, :right2]))
            (_, ints, ligs) = edgestats[(a, b)]
            statsmatrix[i, 1:hist_bins] .= histo(ints, mi1, ma1, hist_bins, true)
            statsmatrix[i, (hist_bins+1):(2*hist_bins)] .= histo(ligs, mi1, ma1, hist_bins, true)
            statsmatrix[i, (2*hist_bins+1):(3*hist_bins)] .= histo(ints, mi2, ma2, hist_bins, false)
            statsmatrix[i, (3*hist_bins+1):(4*hist_bins)] .= histo(ligs, mi2, ma2, hist_bins, false)
        end
        return sort!(stats_df, :nb_ints; rev=true)
    else
        throw(AssertionError("output has to be one of :edges, :nodes, :multi or :stats"))
    end
end