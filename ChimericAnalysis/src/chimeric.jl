struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_selfchimeric=>Int[], :nb_unclassified=>Int[],
        :nb_ints=>Int[], :nb_ints_src=>Int[], :nb_ints_dst=>Int[], :nb_partners=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :meanlen1=>Float64[], :meanlen2=>Float64[],
        :nms1=>Float64[], :nms2=>Float64[])
    edgestats = Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}()
    Interactions(nodes, edges, edgestats, Symbol[])
end

Interactions(alignments::AlignedReads; replicate_id=:first, min_distance=1000, max_ligation_distance=5, filter_types=[], allow_self_chimeras=false) =
    append!(Interactions(), alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                filter_types=filter_types, allow_self_chimeras=allow_self_chimeras)

"""
Load Interactions struct from jld2 file.
"""
Interactions(filepath::String) = load(filepath, "interactions")

Base.length(interactions::Interactions) = nrow(interactions.edges)
function Base.empty!(interactions::Interactions)
    empty!(interactions.nodes)
    empty!(interactions.edges)
    empty!(interactions.edgestats)
    empty!(interactions.replicate_ids)
end

"""
Method of write function which saves the Interactions struct in a jld2 file.
"""
function Base.write(filepath::String, interactions::Interactions)
    if !endswith(filepath, ".jld2")
        throw(ArgumentError("Append '.jld2' to filepath"))
    else
        save(filepath, "interactions", interactions)
    end
end

struct MergedAlignedRead
    alnread::AlignedRead
    pindexpairs::Vector{Tuple{Int,Int}}
    merged::Vector{Bool}
end

MergedAlignedRead(alnread::AlignedRead) = MergedAlignedRead(alnread, Tuple{Int, Int}[], Bool[])

Base.length(mergedread::MergedAlignedRead) = length(mergedread.pindexpairs)

canmerge(alns::AlignedReads, i1::Int, i2::Int) = alns.annotated[i2] && (alns.reads[i1]!==alns.reads[i2]) && (alns.annames[i1]===alns.annames[i2]) &&
((alns.strands[i1] == RNASeqTools.STRAND_POS) ? (alns.leftpos[i1]<=alns.leftpos[i2]) : (alns.rightpos[i1]>=alns.rightpos[i2]))

function mergeparts!(mergedread::MergedAlignedRead, alnread::AlignedRead)
    alns = alnread.alns
    index = @view alns.pindex[alnread.range]
    length(mergedread.pindexpairs) >= length(index) || resize!(mergedread.pindexpairs, length(index))
    length(mergedread.merged) >= length(index) || resize!(mergedread.merged, length(index))
    mergedread.merged .= false
    c::Int = 0
    for (i::Int, ii::Int) in enumerate(index)
        (mergedread.merged[i] || !alns.annotated[ii]) && continue
        c += 1
        nextolpi = findnext(x->canmerge(alns, ii, x), index, i+1)
        if isnothing(nextolpi)
            mergedread.pindexpairs[c] = (ii, ii)
        else
            mergedread.pindexpairs[c] = (ii, index[nextolpi])
            mergedread.merged[nextolpi] = true
        end
    end
    resize!(mergedread.pindexpairs, c)
    return mergedread
end

function hasligationpoint(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; max_distance=5)
    alns = mergedread.alnread.alns
    return ((alns.reads[first(pair1)]==alns.reads[first(pair2)]) && (alns.read_leftpos[first(pair2)]-alns.read_rightpos[first(pair1)] <= max_distance)) ||
    ((alns.reads[last(pair1)]==alns.reads[last(pair2)]) && (alns.read_leftpos[last(pair2)]-alns.read_rightpos[last(pair1)] <= max_distance))
end

function distance(l1::Int, r1::Int, l2::Int, r2::Int; check_order=false)::Float64
    check_order && l1 > l2 && return Inf
    l2>r1 && return l2-r1+1
    l1>r2 && return l1-r2+1
    return 0
end

function ischimeric(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; min_distance=1000, check_annotation=true, check_order=false)
    check_annotation && (mergedread.alnread.alns.annames[first(pair1)] == mergedread.alnread.alns.annames[first(pair2)]) && return false
    (mergedread.alnread.alns.strands[first(pair1)] != mergedread.alnread.alns.strands[first(pair2)]) && return true
    return distance(
        mergedread.alnread.alns.leftpos[first(pair1)],
        mergedread.alnread.alns.rightpos[last(pair1)],
        mergedread.alnread.alns.leftpos[first(pair2)],
        mergedread.alnread.alns.rightpos[last(pair2)]; check_order=check_order) > min_distance
end

function ischimeric(mergedread::MergedAlignedRead; min_distance=1000, check_annotation=true, check_order=false)
    length(mergedread) > 1 || return false
    for (p1, p2) in combinations(mergedread.pindexpairs, 2)
        ischimeric(mergedread, p1, p2; min_distance=min_distance, check_annotation=check_annotation, check_order=check_order) && return true
    end
    return false
end

isselfchimeric(mergedread::MergedAlignedRead) =
    (length(mergedread) == 1) &&
    (mergedread.alnread.alns.leftpos[mergedread.pindexpairs[1][1]] > mergedread.alnread.alns.rightpos[mergedread.pindexpairs[1][2]])

issingle(mergedread::MergedAlignedRead) = (length(mergedread) == 1) && (mergedread.pindexpairs[1][1] != mergedread.pindexpairs[1][2])

function myhash(alns::AlignedReads, i::Int; use_type=true)
    return use_type ? hash(alns.annames[i], hash(alns.antypes[i])) : hash(alns.annames[i])
end

reflen(alns::AlignedReads, pindexpair::Tuple{Int,Int}) =
    max(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)]) -
    min(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)])

function Base.append!(interactions::Interactions, alignments::AlignedReads, replicate_id::Symbol;
                        min_distance=1000, max_ligation_distance=5, filter_types=[], allow_self_chimeras=true)
    if !(String(replicate_id) in interactions.replicate_ids)
        interactions.edges[:, replicate_id] = zeros(Int, nrow(interactions.edges))
        push!(interactions.replicate_ids, replicate_id)
    end
    trans = Dict{UInt, Int}(interactions.nodes[i, :hash]=>i for i in 1:nrow(interactions.nodes))
    edgestats = interactions.edgestats
    exclude_count = 0
    single_count = 0
    chimeric_count = 0
    multi_count = 0
    total_count = 0
    mergedread = MergedAlignedRead(first(alignments))
    for alignment in alignments
        total_count += 1
        !isempty(filter_types) && typein(alignment, filter_types) && (exclude_count += 1; continue)
        mergeparts!(mergedread, alignment)
        is_chimeric = ischimeric(mergedread; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras)
        is_chimeric ? chimeric_count += 1 : single_count +=1

        is_multi = is_chimeric ? length(mergedread)>2 : false
        multi_count += is_multi

        for (i1,_) in mergedread.pindexpairs
            h = myhash(alignments, i1)
            if !(h in keys(trans))
                trans[h] = length(trans) + 1
                push!(interactions.nodes, (alignments.annames[i1], alignments.antypes[i1], alignments.refnames[i1], 0, 0, 0, 0, 0, 0, 0, alignments.strands[i1], h))
            end

            if isselfchimeric(mergedread)
                (interactions.nodes[trans[h], :nb_selfchimeric] += 1)
            elseif !is_chimeric
                issingle(mergedread) ? (interactions.nodes[trans[h], :nb_single] += 1) : (interactions.nodes[trans[h], :nb_unclassified] += 1)
            end
        end

        for (pair1, pair2) in combinations(mergedread.pindexpairs, 2)
            ischimeric(mergedread, pair1, pair2; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras) || continue
            a, b = trans[myhash(alignments, first(pair1))], trans[myhash(alignments, first(pair2))]
            interactions.nodes[a, :nb_ints] += 1
            interactions.nodes[a, :nb_ints_src] += 1
            interactions.nodes[b, :nb_ints] += 1
            interactions.nodes[b, :nb_ints_dst] += 1
            if !((a,b) in keys(edgestats))
                edgestats[(a,b)] = (length(edgestats)+1, Dict{Int,Int}(), Dict{Int,Int}(), Dict{Int,Int}(), Dict{Int,Int}())
                push!(interactions.edges, (a, b, 0, 0, 0.0, 0.0, 0.0, 0.0, (0 for i in 1:length(interactions.replicate_ids))...))
            end
            (iindex, leftintcounter, leftligationcounter, rightintcounter, rightligationcounter) = edgestats[(a, b)]
            interactions.edges[iindex, :nb_ints] += 1
            interactions.edges[iindex, :nb_multi] += is_multi
            interactions.edges[iindex, replicate_id] += 1
            leftpos = alignments.strands[last(pair1)] === STRAND_NEG ? alignments.leftpos[last(pair1)] : alignments.rightpos[last(pair1)]
            rightpos = alignments.strands[first(pair2)] === STRAND_NEG ? alignments.rightpos[first(pair2)] : alignments.leftpos[first(pair2)]
            leftcounter, rightcounter = hasligationpoint(mergedread, pair1, pair2; max_distance=max_ligation_distance) ?
                                            (leftligationcounter, rightligationcounter) : (leftintcounter, rightintcounter)
            leftpos in keys(leftcounter) ? (leftcounter[leftpos]+=1) : (leftcounter[leftpos]=1)
            rightpos in keys(rightcounter) ? (rightcounter[rightpos]+=1) : (rightcounter[rightpos]=1)
            for (s,v) in zip((:meanlen1, :meanlen2, :nms1, :nms2),
                            (reflen(alignments, pair1), reflen(alignments, pair2),
                            max(alignments.nms[first(pair1)], alignments.nms[last(pair1)]),
                            max(alignments.nms[first(pair2)], alignments.nms[last(pair2)])))
                interactions.edges[iindex, s] += (v - interactions.edges[iindex, s]) / interactions.edges[iindex, :nb_ints]
            end
        end
    end
    @info "Processed $total_count reads, found $single_count singles, $chimeric_count ($multi_count) chimeras and excluded $exclude_count"
    return interactions
end

function addpvalues!(interactions::Interactions, genome::Genome, random_model_ecdf::ECDF; fisher_exact_tail="right", include_read_identity=true,
                        include_singles=true, check_interaction_distances=(50,20), bp_parameters=(4,5,1,5,6,4))

    if include_read_identity
        ints_between = interactions.edges[!, :nb_ints]
        other_source = interactions.nodes[interactions.edges[!, :src], :nb_ints_src] .- ints_between
        other_target = interactions.nodes[interactions.edges[!, :dst], :nb_ints_dst] .- ints_between
    else
        check_dict = Dict((s,d)=>n for (s,d,n) in eachrow(interactions.edges[!, [:src, :dst, :nb_ints]]))
        ints_between = [(d,s) in keys(check_dict) ? check_dict[(d,s)]+ints : ints for (s,d,ints) in eachrow(interactions.edges[!, [:src, :dst, :nb_ints]])]
        other_source = interactions.nodes[interactions.edges[!, :src], :nb_ints] .- ints_between
        other_target = interactions.nodes[interactions.edges[!, :dst], :nb_ints] .- ints_between
    end

    if include_singles
        other_source .+= interactions.nodes[interactions.edges[!, :src], :nb_single]
        other_source .+= interactions.nodes[interactions.edges[!, :src], :nb_selfchimeric]
        other_target .+= interactions.nodes[interactions.edges[!, :dst], :nb_single]
        other_target .+= interactions.nodes[interactions.edges[!, :dst], :nb_selfchimeric]
    end
    total_other = sum(interactions.edges[!, :nb_ints]) .- ints_between .- other_source .- other_target .+
        (include_singles ? sum(interactions.nodes[!, :nb_single] .+ interactions.nodes[!, :nb_selfchimeric]) : 0)

    odds_ratio = (ints_between .* total_other) ./ (other_target .* other_source)

    tests = FisherExactTest.(ints_between, other_target, other_source, total_other)
    pvalues_fisher = pvalue.(tests; tail=Symbol(fisher_exact_tail))

    adjp_fisher = adjust(PValues(pvalues_fisher), BenjaminiHochberg())
    interactions.edges[:, :odds_ratio] = odds_ratio
    interactions.edges[:, :pvalue] = pvalues_fisher
    interactions.edges[:, :fdr] = adjp_fisher

    interactions.edges[:, :pred_pvalue] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_fdr] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_score] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_cl1] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_cr1] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_cl2] = fill(NaN, nrow(interactions.edges))
    interactions.edges[:, :pred_cr2] = fill(NaN, nrow(interactions.edges))

    scores = Dict((DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
                    (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
                    (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3])
    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
        gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

    max_dist = maximum(abs.(check_interaction_distances))
    for (c, edge_row) in enumerate(eachrow(interactions.edges))
        if !(isnan(edge_row.modelig1) || isnan(edge_row.modelig2))
            i1, i2 = Int(edge_row.modelig1), Int(edge_row.modelig2)
            strand1, strand2 = interactions.nodes[edge_row[:src], :strand], interactions.nodes[edge_row[:dst], :strand]
            ref1, ref2 = interactions.nodes[edge_row[:src], :ref], interactions.nodes[edge_row[:dst], :ref]
            (((i1 + max_dist) > length(genome.chroms[ref1])) || ((i2 + max_dist) > length(genome.chroms[ref2])) ||
                ((i1 - max_dist) < 1) || ((i2 - max_dist) < 1)) && continue

            s1 = strand1=='+' ?
                genome[ref1][(i1-check_interaction_distances[1]+1):(i1-check_interaction_distances[2])] :
                reverse_complement(genome[ref1][(i1+check_interaction_distances[2]):(i1+check_interaction_distances[1]-1)])
            s2 = strand2=='-' ?
                complement(genome[ref2][(i2-check_interaction_distances[1]+1):(i2-check_interaction_distances[2])]) :
                reverse(genome[ref2][(i2+check_interaction_distances[2]):(i2+check_interaction_distances[1]-1)])

            paln = pairalign(LocalAlignment(), s1, s2, model)
            sco = BioAlignments.score(paln)
            interactions.edges.pred_cl1[c] = paln.aln.a.aln.anchors[1].seqpos
            interactions.edges.pred_cr1[c] = paln.aln.a.aln.anchors[end].seqpos
            interactions.edges.pred_cl2[c] = paln.aln.a.aln.anchors[1].refpos
            interactions.edges.pred_cr2[c] = paln.aln.a.aln.anchors[end].refpos
            edge_row.pred_pvalue, edge_row.pred_score = 1-random_model_ecdf(sco), sco
        end
    end

    nan_index = (!).(isnan.(interactions.edges.pred_pvalue))
    interactions.edges.pred_fdr[nan_index] = adjust(PValues(interactions.edges.pred_pvalue[nan_index]), BenjaminiHochberg())
    return interactions
end

cdsposition(feature::Interval{Annotation}) = hasparam(feature, "cds") ? param(feature, "cds", Int) : 0
function addpositions!(interactions::Interactions, features::Features)
    tus = Dict(hash(name(feature), hash(type(feature)))=>(leftposition(feature), rightposition(feature), cdsposition(feature)) for feature in features)
    for (col_featurepos, col_modes, col_rels) in zip([:left1, :right1, :left2, :right2], [:modeint1, :modeint2, :modelig1, :modelig2], [:rel_int1, :rel_int2, :rel_lig1, :rel_lig2])
        interactions.edges[:, col_featurepos] = Vector{Int}(undef, length(interactions))
        interactions.edges[:, col_modes] = Vector{Float64}(undef, length(interactions))
        interactions.edges[:, col_rels] = Vector{Float64}(undef, length(interactions))
    end
    for edge_row in eachrow(interactions.edges)
        (feature1_left, feature1_right) = tus[interactions.nodes[edge_row[:src], :hash]]
        (feature2_left, feature2_right) = tus[interactions.nodes[edge_row[:dst], :hash]]
        isnegative1 = interactions.nodes[edge_row[:src], :strand] === '-'
        isnegative2 = interactions.nodes[edge_row[:dst], :strand] === '-'
        stats = interactions.edgestats[(edge_row[:src], edge_row[:dst])]
        modeint1 = length(stats[2]) > 0 ? argmax(stats[2]) : NaN64
        modelig1 = length(stats[3]) > 0 ? argmax(stats[3]) : NaN64
        modeint2 = length(stats[4]) > 0 ? argmax(stats[4]) : NaN64
        modelig2 = length(stats[5]) > 0 ? argmax(stats[5]) : NaN64
        (rel_int1, rel_lig1) = ((modeint1, modelig1) .- feature1_left) ./ (feature1_right - feature1_left)
        (rel_int2, rel_lig2) = ((modeint2, modelig2) .- feature2_left) ./ (feature2_right - feature2_left)
        isnegative1 && ((rel_int1, rel_lig1) = (1-rel_int1, 1-rel_lig1))
        isnegative2 && ((rel_int2, rel_lig2) = (1-rel_int2, 1-rel_lig2))
        edge_row[[:left1, :right1, :left2, :right2, :modeint1, :rel_int1, :modeint2, :rel_int2, :modelig1, :rel_lig1, :modelig2, :rel_lig2]] =
            round.((feature1_left, feature1_right, feature2_left, feature2_right, modeint1, rel_int1, modeint2, rel_int2, modelig1, rel_lig1, modelig2, rel_lig2); digits=4)
    end
    interactions.nodes[:, :left] = Vector{Int}(undef, nrow(interactions.nodes))
    interactions.nodes[:, :right] = Vector{Int}(undef, nrow(interactions.nodes))
    interactions.nodes[:, :cds] = Vector{Int}(undef, nrow(interactions.nodes))
    for nodes_row in eachrow(interactions.nodes)
        nodes_row[[:left, :right, :cds]] = tus[nodes_row[:hash]]
    end
    return interactions
end

function histo(ints::Dict{Int,Int}, mi::Int, ma::Int, nbins::Int)
    h = zeros(Int, nbins)
    dbin = (ma-mi+1)/nbins
    for i in mi:ma
        if i in keys(ints)
            h[Int(floor((i-mi)/dbin)) + 1] += ints[i]
        end
    end
    return h
end
function asdataframe(interactions::Interactions; output=:edges, min_reads=5, max_fdr=0.05, max_bp_fdr=0.05, hist_bins=100)
    filter_index = (interactions.edges[!, :nb_ints] .>= min_reads) .& (interactions.edges[!, :fdr] .<= max_fdr) .&
                        ((interactions.edges[!, :pred_fdr] .<= max_bp_fdr) .| isnan.(interactions.edges.pred_fdr))
    out_df = interactions.edges[filter_index, :]
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
        out_df[:, :in_libs] = sum(eachcol(out_df[!, interactions.replicate_ids] .!= 0))
        out_columns = [:name1, :type1, :ref1, :strand1,:left1, :right1, :name2, :type2, :ref2, :strand2, :left2, :right2, :nb_ints, :nb_multi, :in_libs, :pvalue, :fdr, :odds_ratio,
        :pred_pvalue, :pred_fdr, :pred_score, :modeint1, :rel_int1, :modelig1, :rel_lig1, :meanlen1, :nms1, :modeint2, :rel_int2, :modelig2, :rel_lig2, :meanlen2, :nms2]
        return sort!(out_df[!, out_columns], :nb_ints; rev=true)
    elseif output === :nodes
        out_nodes = copy(interactions.nodes)
        for (i,row) in enumerate(eachrow(out_nodes))
            row[:nb_ints] = sum(out_df[(out_df.src .== i) .| (out_df.dst .== i), :nb_ints])
            partners = union!(Set(out_df[out_df.src .== i, :dst]), Set(out_df[out_df.dst .== i, :src]))
            row[:nb_partners] = length(partners)
        end
        return sort!(out_nodes[!, [:name, :type, :ref, :nb_single, :nb_selfchimeric, :nb_unclassified, :nb_ints, :nb_partners]], :nb_single; rev=true)
    elseif output === :stats
        edgestats = interactions.edgestats
        statsmatrix = zeros(nrow(out_df), 4*hist_bins)
        for (i, (a,b, mi1, ma1, mi2, ma2)) in enumerate(zip(out_df[!,:src], out_df[!,:dst], out_df[!, :left1], out_df[!, :right1], out_df[!, :left2], out_df[!, :right2]))
            (_, ints1, ligs1, ints2, ligs2) = edgestats[(a, b)]
            statsmatrix[i, 1:hist_bins] .= histo(ints1, mi1, ma1, hist_bins)
            statsmatrix[i, (hist_bins+1):(2*hist_bins)] .= histo(ligs1, mi1, ma1, hist_bins)
            statsmatrix[i, (2*hist_bins+1):(3*hist_bins)] .= histo(ints2, mi2, ma2, hist_bins)
            statsmatrix[i, (3*hist_bins+1):(4*hist_bins)] .= histo(ligs2, mi2, ma2, hist_bins)
        end
        stats_df = DataFrame(statsmatrix, [
            ["$(i)_ints1" for i in 1:hist_bins]...,
            ["$(i)_ints2" for i in 1:hist_bins]...,
            ["$(i)_lig1" for i in 1:hist_bins]...,
            ["$(i)_lig2" for i in 1:hist_bins]...
            ])
        stats_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        stats_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        stats_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        stats_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        stats_df[:, :left1] = out_df[!, :left1]
        stats_df[:, :right1] = out_df[!, :right1]
        stats_df[:, :left2] = out_df[!, :left2]
        stats_df[:, :right2] = out_df[!, :right2]
        stats_df[:, :nb_ints] = out_df[:, :nb_ints]
        return sort!(stats_df, :nb_ints; rev=true)
    else
        throw(AssertionError("output has to be one of :edges, :nodes, :stats"))
    end
end

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

function chimeric_analysis(features::Features, bams::SingleTypeFiles, results_path::String, conditions::Dict{String, Vector{Int}}, genome::Genome;
                            filter_types=["rRNA", "tRNA"], min_distance=1000, prioritize_type="sRNA", min_prioritize_overlap=0.8, max_bp_fdr=0.05,
                            overwrite_type="IGR", max_ligation_distance=5, is_reverse_complement=true, check_interaction_distances=(45,-10),
                            include_secondary_alignments=true, include_alternative_alignments=false, min_reads=5, max_fdr=0.05, fisher_exact_tail="right",
                            overwrite_existing=false, include_read_identity=true, include_singles=true, allow_self_chimeras=true, position_distribution_bins=50,
                            bp_parameters=(4,5,1,5,6,4), n_genome_samples=200000, fisher_fdr_levels=[0.05,0.1,0.2], bp_fdr_levels=[0.1,0.2,0.5])

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

        seq_length = sum(abs.(check_interaction_distances))

        genome_model_ecdf = ecdf([BioAlignments.score(pairalign(LocalAlignment(),
            i1 % 2 == 0 ? genome.seq[i1:i1+seq_length] : reverse_complement(genome.seq[i1:i1+seq_length]),
            i2 % 2 == 0 ? reverse(genome.seq[i2:i2+seq_length]) : complement(genome.seq[i2:i2+seq_length]),
            model; score_only=true)) for (i1, i2) in eachrow(rand(1:(length(genome.seq)-seq_length), (n_genome_samples,2)))])

        randseq_model_ecdf = ecdf([BioAlignments.score(pairalign(LocalAlignment(),
            randdnaseq(seq_length),
            randdnaseq(seq_length),
            model; score_only=true)) for i in 1:n_genome_samples])

        @info "Using $(summarize(features))"
        isdir(joinpath(results_path, "tables")) || mkpath(joinpath(results_path, "tables"))
        isdir(joinpath(results_path, "jld")) || mkpath(joinpath(results_path, "jld"))
        isdir(joinpath(results_path, "plots")) || mkpath(joinpath(results_path, "plots"))

        interactions = Interactions()
        for (condition, r) in conditions
            @info "Collecting $(length(r)) samples for condition $condition:"
            if (!overwrite_existing && isfile(joinpath(results_path, "tables", "genes_$(condition).csv")) &&
                isfile(joinpath(results_path, "tables", "interactions_$(condition).csv")) &&
                isfile(joinpath(results_path, "tables", "ligation_points_$(condition).csv")) &&
                isfile(joinpath(results_path, "jld", "$(condition).jld2")))

                @info "Found results files. Skipping condition $condition..."
                continue
            end
            replicate_ids = Vector{Symbol}()
            empty!(interactions)
            interactions = Interactions()
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
                    filter_types=filter_types, allow_self_chimeras=allow_self_chimeras)
                empty!(alignments)
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

            p = bp_score_dist_plot(interactions, genome_model_ecdf, randseq_model_ecdf, 0.1)
            save(joinpath(results_path, "plots", "$(condition)_bp_scores_dist.png"), p)

            p2, p3 = bp_clipping_dist_plots(interactions, check_interaction_distances, fisher_fdr_levels, bp_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_clippings_bp.png"), p2)
            save(joinpath(results_path, "plots", "$(condition)_clippings_fisher.png"), p3)

            p4 = log_chimeric_counts_distribution_plot(interactions, fisher_fdr_levels, bp_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_count_dist.png"), p4)

            p5 = log_odds_ratio_distribution_plot(interactions, fisher_fdr_levels, bp_fdr_levels)
            save(joinpath(results_path, "plots", "$(condition)_odds_ratio_dist.png"), p5)
        end
        if !(!overwrite_existing && isfile(joinpath(results_path, "singles.xlsx")) && isfile(joinpath(results_path, "interactions.xlsx")))
            @info "Finalizing..."
            singles = CsvFiles(joinpath(results_path, "tables"); prefix="genes")
            ints = CsvFiles(joinpath(results_path, "tables"); prefix="interactions")
            write(joinpath(results_path, "genes.xlsx"), singles)
            write(joinpath(results_path, "interactions.xlsx"), ints)
        end
        @info "Done."
    end
end