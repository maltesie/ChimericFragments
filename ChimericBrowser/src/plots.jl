function bp_score_dist_plot(interact::Interactions, genome_model_ecdf::ECDF, randseq_model_ecdf::ECDF, plot_fdr_level::Float64)

    interact_ecdf = ecdf([t[6] for t in values(interact.bpstats)])

    max_score = Int(floor(max(
        maximum(interact_ecdf.sorted_values),
        maximum(randseq_model_ecdf.sorted_values),
        maximum(genome_model_ecdf.sorted_values)
    )))
    min_score = Int(floor(min(
        minimum(interact_ecdf.sorted_values),
        minimum(randseq_model_ecdf.sorted_values),
        minimum(genome_model_ecdf.sorted_values)
    )))

    nan_index = .!isnan.(interact.edges.bp_fdr)
    si_fdr = sort(collect(enumerate(zip(interact.edges.bp_fdr[nan_index], interact.edges.bp_pvalue[nan_index]))), by=x->x[2][2], rev=true)

    genome_model_pdf = diff(genome_model_ecdf.(min_score:(max_score+1)))
    randseq_model_pdf = diff(randseq_model_ecdf.(min_score:(max_score+1)))
    interact_pdf = diff(interact_ecdf.(min_score:(max_score+1)))

    p_random = scatter(x=min_score:max_score, y=randseq_model_pdf, fill="tozeroy", name="random")
    p_random_genome = scatter(x=min_score:max_score, y=genome_model_pdf, fill="tozeroy", name="random, from genome")
    p_interact = scatter(x=min_score:max_score, y=interact_pdf, fill="tozeroy", name="around ligation points")
    p = plot([p_random, p_random_genome, p_interact], Layout(yaxis_title="empirical density", xaxis_title="basepairing predictions score distribution"))

    ymax = max(maximum(randseq_model_pdf), maximum(genome_model_pdf), maximum(interact_pdf))
    fdr_p_index = findfirst(x->x[2][1]<=plot_fdr_level, si_fdr)
    fdr_p = isnothing(fdr_p_index) ? 1.0 : si_fdr[fdr_p_index][2][2]
    idx = findfirst(x->((1-genome_model_ecdf(x))<fdr_p), min_score:(max_score+1))
    genome_model_fdr_score = (isnothing(idx) ? 0 : idx - 1) + min_score
    add_shape!(p, line(x0=genome_model_fdr_score, y0=0, x1=genome_model_fdr_score, y1=ymax, line=attr(color="RoyalBlue", width=3), name="FDR=$plot_fdr_level"))

    return p
end

function bp_clipping_dist_plots(interact::Interactions, plotting_fdr_level::Float64)

    l1 = [t[2] for t in values(interact.bpstats)]
    r1 = [t[3] for t in values(interact.bpstats)]
    l2 = [t[4] for t in values(interact.bpstats)]
    r2 = [t[5] for t in values(interact.bpstats)]

    bp_fdrs = adjust(PValues([t[1] for t in values(interact.bpstats)]), BenjaminiHochberg())
    bp_index = bp_fdrs .<= plotting_fdr_level

    h1 = histogram(x=l1[.!bp_index], name="RNA1, left")
    h2 = histogram(x=r1[.!bp_index], name="RNA1, right")
    h3 = histogram(x=l2[.!bp_index], name="RNA2, left")
    h4 = histogram(x=r2[.!bp_index], name="RNA2, right")

    ah1 = histogram(x=l1[bp_index], name="RNA1, left")
    ah2 = histogram(x=r1[bp_index], name="RNA1, right")
    ah3 = histogram(x=l2[bp_index], name="RNA2, left")
    ah4 = histogram(x=r2[bp_index], name="RNA2, right")

    return plot([ah1, ah3, ah2, ah4], Layout(title="alignment clipping for fdr <= $plotting_fdr_level", xaxis_title="position in alignment", yaxis_title="count")),
            plot([h1, h3, h2, h4], Layout(title="alignment clipping for fdr > $plotting_fdr_level", xaxis_title="position in alignment", yaxis_title="count"))
end

function odds_dist_plots(interact::Interactions, plotting_fdr_level::Float64)

    valid_index = interact.edges.odds_ratio .>= 0.0
    logodds = log.(interact.edges.odds_ratio[valid_index])
    infindex = isfinite.(logodds)
    h1 = histogram(x=logodds[infindex], name="all")
    sigfish_index = interact.edges.fisher_fdr[valid_index] .<= plotting_fdr_level
    h2 = histogram(x=logodds[sigfish_index .& infindex], name="fdr <= $plotting_fdr_level")

    all_index = interact.edges.bp_fdr[valid_index] .<= 1.0
    h3 = histogram(x=logodds[all_index], name="all")
    sigpred_index = interact.edges.bp_fdr[valid_index] .<= plotting_fdr_level
    h4 = histogram(x=logodds[sigpred_index .& infindex], name="fdr <= $plotting_fdr_level")

    return plot([h1, h2], Layout(title="Odds ratio distribution for fisher exact test", xaxis_title="log odds ratio", yaxis_title="count")),
            plot([h3, h4], Layout(title="Odds ratio distribution for bp prediction test", xaxis_title="log odds ratio", yaxis_title="count"))
end

function degrees(interact::Interactions; index=trues(nrow(interact.edges)))
    src, dst = interact.edges.src[index], interact.edges.dst[index]
    [length(union!(Set(dst[src .== i]), Set(src[dst .== i]))) for i in unique(vcat(src, dst))]
end

function degree_dist_plots(interact::Interactions, plotting_fdr_level::Float64)

    sigfish_index = interact.edges.fisher_fdr .<= plotting_fdr_level
    degs = counter(degrees(interact; index=sigfish_index))
    h2 = scatter(x=log.(collect(keys(degs))), y=log.(collect(values(degs)./sum(values(degs)))),
        name="fisher fdr <= $plotting_fdr_level", mode="markers")

    sigpred_index = interact.edges.bp_fdr .<= plotting_fdr_level
    degs = counter(degrees(interact; index=sigpred_index))
    h4 = scatter(x=log.(collect(keys(degs))), y=log.(collect(values(degs)./sum(values(degs)))),
        name="bp fdr <= $plotting_fdr_level", mode="markers")

    plot(h2, Layout(title="Node degree dristribution for fisher fdr <= $plotting_fdr_level", xaxis_title="log degree", yaxis_title="log ratio")),
        plot(h4, Layout(title="Node degree dristribution for bp fdr <= $plotting_fdr_level", xaxis_title="log degree", yaxis_title="log ratio"))
end

function annotation_type_heatmap(interact::Interactions, plotting_fdr_level::Float64)

    type1 = interact.nodes.type[interact.edges.src]
    type2 = interact.nodes.type[interact.edges.dst]
    types = collect(union!(Set(type1), Set(type2)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))

    types_counter = zeros(Float64, (length(types), length(types)))
    type1 = interact.nodes.type[interact.edges.src[interact.edges.fisher_fdr .<= plotting_fdr_level]]
    type2 = interact.nodes.type[interact.edges.dst[interact.edges.fisher_fdr .<= plotting_fdr_level]]
    for (t1, t2) in zip(type1, type2)
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    cross_types = ["$t1-$t2" for (i,t1) in enumerate(types) for t2 in types[i:end]]
    cross_values =  [types_counter[type_trans[t1], type_trans[t2]] for (i,t1) in enumerate(types) for t2 in types[i:end] ]
    b1 = bar(x=cross_types, y=cross_values, name="fisher fdr <= $plotting_fdr_level")

    types_counter = zeros(Float64, (length(types), length(types)))
    type1 = interact.nodes.type[interact.edges.src[interact.edges.bp_fdr .<= plotting_fdr_level]]
    type2 = interact.nodes.type[interact.edges.dst[interact.edges.bp_fdr .<= plotting_fdr_level]]
    for (t1, t2) in zip(type1, type2)
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    cross_values =  [types_counter[type_trans[t1], type_trans[t2]] for (i,t1) in enumerate(types) for t2 in types[i:end] ]
    b2 = bar(x=cross_types, y=cross_values, name="bp fdr <= $plotting_fdr_level")

    singles = [sum(interact.nodes.nb_single[interact.nodes.type .== t]) for t in types]
    b3 = bar(x=types, y=singles)

    return plot([b1, b2], Layout(title="Interactions between annotation types", yaxis_title="count")),
            plot(b3, Layout(title="Single reads per annotation type", yaxis_title="count"))
end

function plot_pair(interact::Interactions, t::String, max_fdr::Float64)

    p1, p2 = if t == "bp"
        bp_clipping_dist_plots(interact, max_fdr)
    elseif t == "odds"
        odds_dist_plots(interact, max_fdr)
    elseif t == "annotation"
        annotation_type_heatmap(interact, max_fdr)
    elseif t == "degree"
        degree_dist_plots(interact, max_fdr)
    else
        bp_clipping_dist_plots(interact, max_fdr)
    end
    return p1, p2
end

alnchar(x::DNA, y::DNA) =
    if (x == DNA_A && y == DNA_T) || (x == DNA_T && y == DNA_A) || (x == DNA_C && y == DNA_G) || (x == DNA_G && y == DNA_C)
        '|'
    elseif (x == DNA_G && y == DNA_T) || (x == DNA_T && y == DNA_G)
        '·'
    else
        ' '
    end
function basepairing_string(aln::PairwiseAlignment, n1::String, n2::String, offset1::Int, offset2::Int, al1::Int, ar1::Int, al2::Int, ar2::Int)
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    posw = max(ndigits(offset1 + anchors[end].seqpos)+1, ndigits(offset2 - anchors[1].refpos)+1, ndigits(al1), ndigits(al2)) + 1

    seqpos = offset1 + anchors[1].seqpos + 1
    refpos = offset2 - anchors[1].refpos + 2
    seqbuf = IOBuffer()
    refbuf = IOBuffer()
    head1buf = IOBuffer()
    head2buf = IOBuffer()
    matbuf = IOBuffer()

    print(seqbuf, lpad(seqpos>0 ? "+$seqpos " : "$(seqpos-1) ", posw))
    print(refbuf, lpad(refpos>0 ? "+$refpos " : "$(refpos-1) ", posw))
    print(head1buf, lpad("$al1 ", posw))
    print(head2buf, lpad("$al2 ", posw))
    #print(matbuf, " "^posw)

    next_xy = iterate(aln)
    while next_xy !== nothing
        (x, y), s = next_xy
        next_xy = iterate(aln ,s)

        if x != gap(eltype(seq))
            seqpos += 1
        end
        if y != gap(eltype(ref))
            refpos -= 1
        end

        print(seqbuf, RNA(x))
        print(refbuf, RNA(y))
        print(matbuf, alnchar(x, y))
    end

    seqpos -= 1
    refpos += 1
    print(seqbuf, seqpos > 0 ? " +$seqpos" : " $(seqpos-1)")
    print(refbuf, refpos > 0 ? " +$refpos" : " $(refpos-1)")

    seqs = String(take!(seqbuf))
    mats = String(take!(matbuf))
    refs = String(take!(refbuf))

    noffset = Int(floor(length(mats)/2))
    print(head1buf, rpad(lpad(n1, length(n1)+noffset-Int(floor(length(n1)/2))), length(mats)))
    print(head2buf, rpad(lpad(n2, length(n2)+noffset-Int(floor(length(n2)/2))), length(mats)))
    print(head1buf, " $ar1")
    print(head2buf, " $ar2")
    head1s = String(take!(head1buf))
    head2s = String(take!(head2buf))

    head1s * "<br>" * seqs * "<br>" * lpad(mats, length(mats) + posw)  * "<br>" * refs * "<br>" * head2s
end
function alignment_ascii_plot(i1::Int, i2::Int, p1::Int, p2::Int, interact::Interactions,
        genome::Dict{String, BioSequences.LongDNA{4}}, check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel)

    n1::String, ref1::String, strand1::Char, l1::Int, r1::Int, c1::Int = interact.nodes[i1, [:name, :ref, :strand, :left, :right, :cds]]
    n2::String, ref2::String, strand2::Char, l2::Int, r2::Int, c2::Int = interact.nodes[i2, [:name, :ref, :strand, :left, :right, :cds]]

    _, al1, ar1, al2, ar2, _ = interact.bpstats[(i1, p1, i2, p2)]

    al1 = p1 + (strand1=='-' ? 1 : -1) * (check_interaction_distances[1] - al1)
    ar1 = p1 + (strand1=='-' ? 1 : -1) * (check_interaction_distances[1] - ar1)
    al2 = p2 + (strand2=='-' ? -1 : 1) * (check_interaction_distances[1] - al2)
    ar2 = p2 + (strand2=='-' ? -1 : 1) * (check_interaction_distances[1] - ar2) - 1

    s1 = strand1=='+' ?
        genome[ref1][(p1-check_interaction_distances[1]):(p1-check_interaction_distances[2])] :
        BioSequences.reverse_complement(genome[ref1][(p1+check_interaction_distances[2]):(p1+check_interaction_distances[1])])
    s2 = strand2=='-' ?
        BioSequences.complement(genome[ref2][(p2-check_interaction_distances[1]):(p2-check_interaction_distances[2])]) :
        BioSequences.reverse(genome[ref2][(p2+check_interaction_distances[2]):(p2+check_interaction_distances[1])])
    p = pairalign(LocalAlignment(), s1, s2, model)

    basepairing_string(alignment(p), n1, n2,
        strand1=='+' ? ((p1-check_interaction_distances[1])-(c1>0 ? c1 : l1)) : ((c1>0 ? c1 : r1)-(p1+check_interaction_distances[1])),
        strand2=='-' ? ((c2>0 ? c2 : r2)-(p2-check_interaction_distances[1])) : ((p2+check_interaction_distances[1])-(c2>0 ? c2 : l2)),
        al1, ar1, al2, ar2)
end

normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=0, to_max=100) = Int(floor(to_min + value * (to_max-to_min)))
function cdsframestring(p::Int, idx::Int, interact::Interactions)
    cds, left, right = interact.nodes[idx, [:cds, :left, :right]]
    tp = interact.nodes.strand[idx] == '-' ? (cds > 0 ? cds : right)-p+1 : p-(cds > 0 ? cds : left)+1
    tp <= 0 && (tp -= 1)
    return tp > 0 ? "+$tp" : "$tp"
end

function edge_figure(edge_data::Dash.JSON3.Object, interact::Interactions, genome::Dict{String,BioSequences.LongDNA{4}},
        check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel, max_fdr::Float64)
    src, dst = parse(Int, edge_data["source"]), parse(Int, edge_data["target"])
    name1, name2 = interact.nodes[src, :name], interact.nodes[dst, :name]
    ((src,dst) in keys(interact.edgestats)) || return no_ligs_edge_figure(edge_data["interactions"])
    fdrs = adjust(PValues([interact.bpstats[(src, i1, dst, i2)][1] for (i1, i2) in keys(interact.edgestats[(src,dst)][3])]), BenjaminiHochberg())
    fdr_keys = [p for (i, p) in enumerate(keys(interact.edgestats[(src,dst)][3])) if fdrs[i] <= max_fdr]
    isempty(fdr_keys) && return no_ligs_edge_figure(edge_data["interactions"])
    points1, points2 = first.(fdr_keys), last.(fdr_keys)
    tickpos1, tickpos2 = [minimum(points1), maximum(points1)], [minimum(points2), maximum(points2)]
    ticks1, ticks2 = [cdsframestring(p, src, interact) for p in tickpos1], [cdsframestring(p, dst, interact) for p in tickpos2]
    maxints = maximum(values(interact.edgestats[(src, dst)][3]))
    sizes = [ceil(interact.edgestats[(src, dst)][3][p]/maxints*4)*5 for p in zip(points1, points2)]
    colors = [fdr for fdr in fdrs if fdr <= max_fdr]
    bp_plots = [alignment_ascii_plot(src,dst,p1,p2,interact,genome,check_interaction_distances, model) *
        "<br><br>ligation point: ($(cdsframestring(p1, src, interact)), $(cdsframestring(p2, dst, interact)))<br>FDR: $(round(c, digits=4))<br>" *
        "supporting reads: $(interact.edgestats[(src,dst)][3][(p1,p2)])" for (p1,p2,c) in zip(points1, points2, colors)]
    return plot(scatter(y=points1, x=points2, mode="markers", marker=attr(size=sizes, color=colors, colorbar=attr(title="FDR", orientation="v"),
            colorscale="Reds", reversescale=true, cmin=0.0, cmax=1.0), name = "ligation points",
            text=bp_plots, hoverinfo="text", hovertemplate="<span style=\"font-family:'Lucida Console', monospace\">%{text}</span><extra></extra>"),
        Layout(title = "$(edge_data["interactions"]) interactions, $(sum(interact.edgestats[(src,dst)][3][p] for p in zip(points1, points2))) ligation points",
            yaxis=attr(title="RNA1: $name1", tickmode="array", tickvals=tickpos1, ticktext=ticks1),
            xaxis=attr(title="RNA2: $name2", tickmode="array", tickvals=tickpos2, ticktext=ticks2)))
end

function count_ligation_sites_as1(nodeid::Int64, node_data::Dash.JSON3.Object, interact::Interactions, bp_len::Int, max_fdr::Float64)
    counts = Dict{Int, Dict{Int,Int}}()
    isnegative = interact.nodes.strand[nodeid] == '-'
    partners = node_data["lig_as_rna1"]
    for partner in partners
        ligation_points = keys(interact.edgestats[(nodeid, partner)][3])
        fdr = adjust(PValues([interact.bpstats[(nodeid, p[1], partner, p[2])][1] for p in ligation_points]), BenjaminiHochberg())
        for ((r1, r2), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            p = (nodeid, r1, partner, r2)
            c = interact.bpstats[p][7]
            for pl in interact.bpstats[p][2]:interact.bpstats[p][3]
                p1 = r1 + (isnegative ? 1 : -1) * (bp_len - pl)
                if p1 in keys(counts)
                    partner in keys(counts[p1]) ? counts[p1][partner] += c : counts[p1][partner] = c
                else
                    counts[p1] = Dict(partner=>c)
                end
            end
        end
    end
    return counts
end
function count_ligation_sites_as2(nodeid::Int64, node_data::Dash.JSON3.Object, interact::Interactions, bp_len::Int, max_fdr::Float64)
    counts = Dict{Int, Dict{Int,Int}}()
    isnegative = interact.nodes.strand[nodeid] == '-'
    partners = node_data["lig_as_rna2"]
    for partner in partners
        ligation_points = keys(interact.edgestats[(partner, nodeid)][3])
        fdr = adjust(PValues([interact.bpstats[(partner, p[1], nodeid, p[2])][1] for p in ligation_points]), BenjaminiHochberg())
        for ((r1, r2), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            p = (partner, r1, nodeid, r2)
            c = interact.bpstats[p][7]
            for pl in interact.bpstats[p][4]:interact.bpstats[p][5]
                p2 = r2 + (isnegative ? -1 : 1) * (bp_len - pl) - 1
                if p2 in keys(counts)
                    partner in keys(counts[p2]) ? counts[p2][partner] += c : counts[p2][partner] = c
                else
                    counts[p2] = Dict(partner=>c)
                end
            end
        end
    end
    return counts
end
function nested_join(s::Vector{String}, n::Int, max_len::Int, it_minor::String, it_major::String)
    mys = length(s) > max_len ? [s[1:min(length(s), max_len)]..., "..."] : s
    join([join(mys[((i-1)*n+1):min((i*n), length(mys))], it_minor) for i in 1:(1+div(length(mys), n))], it_major)
end
function node_figure(node_data::Dash.JSON3.Object, interact::Interactions, interaction_distances::Tuple{Int64, Int64}, max_fdr::Float64)
    idx = parse(Int, node_data["id"])
    name = interact.nodes.name[idx]
    lp_rna1 = count_ligation_sites_as1(idx, node_data, interact, interaction_distances[1], max_fdr)
    lp_rna2 = count_ligation_sites_as2(idx, node_data, interact, interaction_distances[1], max_fdr)
    length(lp_rna1) == 0 && length(lp_rna2) == 0 && return empty_node_figure
    minpos = minimum(minimum(k for k in keys(lp)) for lp in (lp_rna1, lp_rna2) if length(lp)>0)
    maxpos = maximum(maximum(k for k in keys(lp)) for lp in (lp_rna1, lp_rna2) if length(lp)>0)
    tickpos = [minpos, maxpos]
    ticktext = [cdsframestring(p, idx, interact) for p in tickpos]
    return plot([
            begin
                kv = sort(collect(ligationpoints), by=x->x[1])
                positions, counts = first.(kv), sum.(values.(last.(kv)))
                allpositions = length(positions) > 0 ? collect(minimum(positions):maximum(positions)) : Int[]
                pindex = in.(allpositions, Ref(positions))
                indextrans = [pindex[i] ? sum(view(pindex, 1:i)) : 0 for i in eachindex(allpositions)]
                allcounts = [pindex[i] ? counts[indextrans[i]] : 0 for i in eachindex(allpositions)]
                nb_binding = Int[pindex[i] ? length(kv[indextrans[i]][2]) : 0 for i in eachindex(allpositions)]
                partners = [pindex[i] ? nested_join(["$(interact.nodes.name[n]): $c" for (n, c) in sort(collect(ligationpoints[p]),
                    by=x->x[2], rev=true)], 3, 17, ", ", "<br>") : "" for (i, p) in enumerate(allpositions)]
                ticks = [cdsframestring(p, idx, interact) for p in tickpos]
                hover_texts = ["position: $(cdsframestring(p, idx, interact)) ($p)<br># of partners here: $b<br>total reads count: $c<br>$t"
                    for (p, c, t, b) in zip(allpositions, allcounts, partners, nb_binding)]
                scatter(x = allpositions, y = allcounts, fill="tozeroy", name=legend, text=hover_texts, hoverinfo="text")
            end
            for (ligationpoints, legend) in zip((lp_rna1, lp_rna2), ("as RNA1", "as RNA2"))
        ],
        Layout(title = "$name on $(interact.nodes.ref[idx]) ($(interact.nodes.strand[idx]))", xaxis_title = "position", yaxis_title = "count",
            xaxis=attr(tickmode="array", tickvals=tickpos, ticktext=ticktext), showlegend=false)
    )
end

const empty_figure = (
    data = [
        (x = [], y = [], type = "scatter", name = ""),
    ],
    layout = (title = "Please select an edge<br>or a node in the graph.",)
)

const empty_node_figure = (
    data = [
        (x = [], y = [], type = "scatter", name = ""),
    ],
    layout = (title = "No basepairing predictions found<br>for selected FDR.",)
)

no_ligs_edge_figure(nb_ints::Int) = (
    data = [
        (x = [], y = [], type = "scatter", name = ""),
    ],
    layout = (title = "$nb_ints interactions.<br>No ligation points found.",)
)