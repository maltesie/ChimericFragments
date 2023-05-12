function bp_score_dist_plot(interact::Interactions, genome_model_ecdf::ECDF, randseq_model_ecdf::ECDF, plot_fdr_level::Float64)

    interact_ecdf = ecdf([t[6] for t in values(interact.bpstats)])

    max_score = Int(floor(max(
        maximum(interact_ecdf.sorted_values),
        maximum(randseq_model_ecdf.sorted_values),
        maximum(genome_model_ecdf.sorted_values)
    )))

    nan_index = .!isnan.(interact.edges.bp_fdr)
    si_fdr = sort(collect(enumerate(zip(interact.edges.bp_fdr[nan_index], interact.edges.bp_pvalue[nan_index]))), by=x->x[2][2], rev=true)

    genome_model_pdf = diff(genome_model_ecdf.(1:(max_score+1)))
    randseq_model_pdf = diff(randseq_model_ecdf.(1:(max_score+1)))
    interact_pdf = diff(interact_ecdf.(1:(max_score+1)))

    p_random = scatter(x=1:max_score, y=randseq_model_pdf, fill="tozeroy", name="random")
    p_random_genome = scatter(x=1:max_score, y=genome_model_pdf, fill="tozeroy", name="random, from genome")
    p_interact = scatter(x=1:max_score, y=interact_pdf, fill="tozeroy", name="around ligation points")
    p = plot([p_random, p_random_genome, p_interact], Layout(title="basepairing predictions score distribution"))

    ymax = max(maximum(randseq_model_pdf), maximum(genome_model_pdf), maximum(interact_pdf))
    fdr_p_index = findfirst(x->x[2][1]<=plot_fdr_level, si_fdr)
    fdr_p = isnothing(fdr_p_index) ? 1.0 : si_fdr[fdr_p_index][2][2]
    idx = findfirst(x->((1-genome_model_ecdf(x))<fdr_p), 1:(max_score+1))
    genome_model_fdr_score = isnothing(idx) ? 0 : idx - 1
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

    return plot([ah1, ah3, ah2, ah4], Layout(title="alignment clipping for fdr <= $plotting_fdr_level")),
            plot([h1, h3, h2, h4], Layout(title="alignment clipping for fdr > $plotting_fdr_level"))
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

    plot([h1, h2]), plot([h3, h4])
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

    plot(h2), plot(h4)
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
    types_counter ./= sum(types_counter)
    h1 = heatmap(x=types, y=types, z=copy(types_counter))

    types_counter = zeros(Float64, (length(types), length(types)))
    type1 = interact.nodes.type[interact.edges.src[interact.edges.bp_fdr .<= plotting_fdr_level]]
    type2 = interact.nodes.type[interact.edges.dst[interact.edges.bp_fdr .<= plotting_fdr_level]]
    for (t1, t2) in zip(type1, type2)
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    types_counter ./= sum(types_counter)
    h2 = heatmap(x=types, y=types, z=copy(types_counter))

    return plot(h1), plot(h2)
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
        '⋅'
    else
        ' '
    end
function basepairing_string(aln::PairwiseAlignment, n1::String, n2::String, offset1::Int, offset2::Int, al1::Int, ar1::Int, al2::Int, ar2::Int)
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    posw = max(ndigits(offset1 + anchors[end].seqpos)+1, ndigits(offset2 - anchors[1].refpos)+1, ndigits(al1), ndigits(ar1)) + 1

    i = 0
    seqpos = offset1 + anchors[1].seqpos
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

        i += 1
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

    _, al1, ar1, al2, ar2, _ = interact.bpstats[(p1,p2)]

    al1 = strand1=='+' ? p1-check_interaction_distances[1]+al1 : p1+check_interaction_distances[1]-al1
    ar1 = strand1=='+' ? p1-check_interaction_distances[1]+ar1 : p1+check_interaction_distances[1]-ar1
    al2 = strand1=='-' ? p2-check_interaction_distances[1]+al2 : p2+check_interaction_distances[1]-al2
    ar2 = strand1=='-' ? p2-check_interaction_distances[1]+ar2 : p2+check_interaction_distances[1]-ar2

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