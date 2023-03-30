struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_selfchimeric=>Int[], :nb_unclassified=>Int[],
        :nb_ints=>Int[], :nb_ints_src=>Int[], :nb_ints_dst=>Int[], :nb_partners=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :meanlen1=>Float64[], :meanlen2=>Float64[], :nms1=>Float64[], :nms2=>Float64[])
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

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

function Base.append!(interactions::Interactions, alignments::AlignedReads, replicate_id::Symbol;
    min_distance=1000, max_ligation_distance=5, filter_types=["rRNA", "tRNA"], allow_self_chimeras=false, is_paired_end=true)

    if !(String(replicate_id) in interactions.replicate_ids)
        interactions.edges[:, replicate_id] = zeros(Int, nrow(interactions.edges))
        push!(interactions.replicate_ids, replicate_id)
    end

    counts = zeros(Int, 6) #single, unclassified, selfchimeric, exclude, chimeric, multi
    trans = Dict{UInt, Int}(interactions.nodes.hash[i]=>i for i in 1:nrow(interactions.nodes))
    edgestats = interactions.edgestats
    node_ints = Matrix{Int}(interactions.nodes[:, [:nb_single, :nb_unclassified, :nb_selfchimeric, :nb_ints, :nb_ints_src, :nb_ints_dst]])
    node_ids = zeros(Int, length(interactions.nodes.hash))
    node_hashs = Vector{UInt}(interactions.nodes.hash)
    nb_nodes_start = nrow(interactions.nodes)
    edge_ints = hcat(Matrix{Int}(interactions.edges[:, [:src, :dst, :nb_ints, :nb_multi]]), zeros(Int, nrow(interactions.edges)))
    edge_floats = Matrix{Float64}(interactions.edges[:, [:meanlen1, :meanlen2, :nms1, :nms2]])
    mergedread = MergedAlignedRead(alignments)
    h = UInt(1)
    for alignment in alignments
        if (!isempty(filter_types) && typein(alignment, filter_types)) || !(hasannotation(alignment))
            counts[4] += 1
            continue
        end
        mergeparts!(mergedread, alignment)
        is_chimeric = ischimeric(mergedread; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras)
        counts[5] += is_chimeric
        is_multi = length(mergedread)>2
        counts[6] += is_multi

        for (i1,_) in mergedread.pindexpairs
            h = myhash(alignments, i1)
            if !(h in keys(trans))
                idx = length(trans) + 1
                trans[h] = idx
                if idx > length(node_hashs)
                    resize!(node_hashs, length(node_hashs) + 100000)
                    resize!(node_ids, length(node_hashs) + 100000)
                    node_ints = vcat(node_ints, zeros(Int, 100000, 6))
                end
                node_ids[idx] = i1
                node_hashs[idx] = h
            end
        end

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

        for (i, pair1) in enumerate(mergedread.pindexpairs), pair2 in (@view mergedread.pindexpairs[i+1:end])
            ischimeric(mergedread, pair1, pair2; min_distance=min_distance, check_annotation=!allow_self_chimeras, check_order=allow_self_chimeras) || continue
            a, b = trans[myhash(alignments, first(pair1))], trans[myhash(alignments, first(pair2))]
            node_ints[a, 4] += 1
            node_ints[a, 5] += 1
            node_ints[b, 4] += 1
            node_ints[b, 6] += 1
            if !((a,b) in keys(edgestats))
                idx = length(edgestats)+1
                edgestats[(a,b)] = (idx, Dict{Int,Int}(), Dict{Int,Int}(), Dict{Int,Int}(), Dict{Int,Int}())
                if idx > size(edge_ints)[1]
                    edge_ints = vcat(edge_ints, zeros(Int, 100000, 5))
                    edge_floats = vcat(edge_floats, zeros(Float64, 100000, 4))
                end
                edge_ints[idx, 1] = a
                edge_ints[idx, 2] = b
            end
            (iindex, leftintcounter, leftligationcounter, rightintcounter, rightligationcounter) = edgestats[(a, b)]
            edge_ints[iindex, 3] += 1
            edge_ints[iindex, 4] += is_multi
            edge_ints[iindex, 5] += 1
            leftpos = alignments.strands[last(pair1)] === STRAND_NEG ? alignments.leftpos[last(pair1)] : alignments.rightpos[last(pair1)]
            rightpos = alignments.strands[first(pair2)] === STRAND_NEG ? alignments.rightpos[first(pair2)] : alignments.leftpos[first(pair2)]
            leftcounter, rightcounter = hasligationpoint(mergedread, pair1, pair2; max_distance=max_ligation_distance) ?
                                            (leftligationcounter, rightligationcounter) : (leftintcounter, rightintcounter)
            leftpos in keys(leftcounter) ? (leftcounter[leftpos]+=1) : (leftcounter[leftpos]=1)
            rightpos in keys(rightcounter) ? (rightcounter[rightpos]+=1) : (rightcounter[rightpos]=1)
            for (i,v) in enumerate((reflen(alignments, pair1), reflen(alignments, pair2),
                            max(alignments.nms[first(pair1)], alignments.nms[last(pair1)]),
                            max(alignments.nms[first(pair2)], alignments.nms[last(pair2)])))
                edge_floats[iindex, i] += (v - edge_floats[iindex, i]) / edge_ints[iindex, 5]
            end
        end
    end

    resize!(interactions.nodes, length(trans))
    nr = (nb_nodes_start+1):length(trans)
    view(interactions.nodes.name, nr) .= view(alignments.annames, view(node_ids, nr))
    view(interactions.nodes.type, nr) .= view(alignments.antypes, view(node_ids, nr))
    view(interactions.nodes.ref, nr) .= view(alignments.refnames, view(node_ids, nr))
    view(interactions.nodes.strand, nr) .= view(alignments.strands, view(node_ids, nr))
    view(interactions.nodes.hash, nr) .= view(node_hashs, nr)

    nr = 1:length(trans)
    for (i, n) in enumerate((:nb_single, :nb_selfchimeric, :nb_unclassified, :nb_ints, :nb_ints_src, :nb_ints_dst))
        interactions.nodes[!, n] .= view(node_ints, nr, i)
    end

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
    infotable = DataFrame(
        :total=>[nread(alignments)], :chimeric=>[counts[5]], :multi=>[counts[6]], :self=>[counts[3]],
        :single=>[counts[1]], :filtered=>[counts[4]], :no_class=>[counts[2]],
    )
    @info "Classification of reads:\n" * DataFrames.pretty_table(String, infotable, nosubheader=true)
    return interactions
end

alnchar(x::DNA, y::DNA) =
    if (x == DNA_A && y == DNA_T) || (x == DNA_T && y == DNA_A) || (x == DNA_C && y == DNA_G) || (x == DNA_G && y == DNA_C)
        '|'
    elseif (x == DNA_G && y == DNA_T) || (x == DNA_T && y == DNA_G)
        '⋅'
    else
        ' '
    end

function basepairing_string(aln::PairwiseAlignment, offset1::Int, offset2::Int; width::Integer=200)
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    # width of position numbers
    posw = ndigits(max(abs(offset1 + anchors[end].seqpos), abs(offset2 - anchors[1].refpos))) + 2
    outstring = ""
    i = 0
    seqpos = anchors[1].seqpos
    refpos = anchors[1].refpos
    seqbuf = IOBuffer()
    refbuf = IOBuffer()
    matbuf = IOBuffer()
    next_xy = iterate(aln)
    while next_xy !== nothing
        (x, y), s = next_xy
        next_xy = iterate(aln ,s)

        i += 1
        if x != gap(eltype(seq))
            seqpos += 1
        end
        if y != gap(eltype(ref))
            refpos += 1
        end

        if i % width == 1
            print(seqbuf, "RNA1:", lpad(seqpos>0 ? "+$seqpos" : "$(seqpos-1)", posw), ' ')
            print(refbuf, "RNA2:", lpad(refpos>0 ? "+$refpos" : "$(refpos-1)", posw), ' ')
            print(matbuf, " "^(posw + 6))
        end

        print(seqbuf, RNA(x))
        print(refbuf, RNA(y))
        print(matbuf, alnchar(x, y))

        if i % width == 0
            print(seqbuf, seqpos > 0 ? " +$seqpos" : " $(seqpos-1)")
            print(refbuf, refpos > 0 ? " +$refpos" : " $(refpos-1)")
            print(matbuf)

            outstring *= String(take!(seqbuf)) * "\n" * String(take!(matbuf)) * "\n" * String(take!(refbuf)) * "\n⋅⋅⋅\n"

            if next_xy !== nothing
                seek(seqbuf, 0)
                seek(matbuf, 0)
                seek(refbuf, 0)
            end
        end
    end

    if i % width != 0
        print(seqbuf, seqpos > 0 ? " +$seqpos" : " $(seqpos-1)")
        print(refbuf, refpos > 0 ? " +$refpos" : " $(refpos-1)")
        print(matbuf)

        outstring *= String(take!(seqbuf)) * "\n" * String(take!(matbuf)) * "\n" * String(take!(refbuf))
    end
    outstring
end

function addpvalues!(interactions::Interactions, genome::Genome, random_model_ecdf::ECDF; fisher_exact_tail="right", include_read_identity=true,
                        include_singles=true, check_interaction_distances=(30,0), bp_parameters=(4,5,0,7,8,3), shift_weight=0.5)

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



    tests = FisherExactTest.(ints_between, other_target, other_source, total_other)
    odds_ratio = [t.ω for t in tests]
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
    complement_genome = Genome(complement(genome.seq), genome.chroms)
    reverse_genome = Genome(copy(genome.seq), genome.chroms)
    for (_, seq) in reverse_genome
        reverse!(seq)
    end
    reverse_complement_genome = Genome(complement(genome.seq), genome.chroms)
    for (_, seq) in reverse_complement_genome
        reverse!(seq)
    end

    for (c, edge_row) in enumerate(eachrow(interactions.edges))
        if !(isnan(edge_row.modelig1) || isnan(edge_row.modelig2))
            i1, i2 = Int(edge_row.modelig1), Int(edge_row.modelig2)
            strand1, strand2 = interactions.nodes[edge_row[:src], :strand], interactions.nodes[edge_row[:dst], :strand]
            ref1, ref2 = interactions.nodes[edge_row[:src], :ref], interactions.nodes[edge_row[:dst], :ref]
            (((i1 + max_dist) > length(genome.chroms[ref1])) || ((i2 + max_dist) > length(genome.chroms[ref2])) ||
                ((i1 - max_dist) < 1) || ((i2 - max_dist) < 1)) && continue

            s1 = if strand1=='+'
                view(genome[ref1], (i1-check_interaction_distances[1]):(i1-check_interaction_distances[2]))
            else
                l = length(genome.chroms[ref1])
                c1, c2 = l+1-(i1+check_interaction_distances[2]), l+1-(i1+check_interaction_distances[1])
                view(reverse_complement_genome[ref1], c2:c1)
            end
            s2 = if strand2=='-'
                view(complement_genome[ref2], (i2-check_interaction_distances[1]):(i2-check_interaction_distances[2]))
            else
                l = length(genome.chroms[ref2])
                c1, c2 = l+1-(i2+check_interaction_distances[2]), l+1-(i2+check_interaction_distances[1])
                view(reverse_genome[ref2], c2:c1)
            end
            paln = pairalign(LocalAlignment(), s1, s2, model)
            sco = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))
            interactions.edges.pred_cl1[c] = paln.aln.a.aln.anchors[1].seqpos + 1
            interactions.edges.pred_cr1[c] = paln.aln.a.aln.anchors[end].seqpos
            interactions.edges.pred_cl2[c] = paln.aln.a.aln.anchors[1].refpos + 1
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