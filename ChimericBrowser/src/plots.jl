hasfdrvalues(ints::Interactions) = "fdr" in names(ints.edges)

function checkinteract(ints::Interactions, verified_pairs::Vector{Tuple{String,String}};
        min_reads=5, max_fdr=0.05, case_sensitive=true)
	verified_dict = merge(Dict(pair=>zeros(5) for pair in verified_pairs),
							Dict(reverse(pair)=>zeros(5) for pair in verified_pairs))
    !case_sensitive && (verified_dict = Dict((uppercase(p[1]), uppercase(p[2]))=>v for (p,v) in verified_dict))
	df = asdataframe(ints; min_reads=min_reads, max_fdr =max_fdr)
    for row in eachrow(df)
        key = (row[:name1], row[:name2])
        !case_sensitive && (key = (uppercase(key[1]), uppercase(key[2])))
        if key in keys(verified_dict)
            verified_dict[key][1] = row[:in_libs]
            verified_dict[key][2] = row[:nb_ints]
            verified_dict[key][3] = row[:fdr]
            verified_dict[key][4] = row[:pred_fdr]
            verified_dict[key][5] = row[:odds_ratio]
        end
    end

	sorted_keys = vcat([[pair, reverse(pair)] for pair in verified_pairs]...)
    !case_sensitive && (sorted_keys = [(uppercase(p[1]), uppercase(p[2])) for p in sorted_keys])
	m = reduce(hcat, [verified_dict[key] for key in sorted_keys])'
    n1s = String[n[1] for n in verified_pairs]
    n2s = String[n[2] for n in verified_pairs]
	verified_stats = DataFrame(
		name1=vec(vcat(n1s', n2s')),
		name2=vec(vcat(n2s', n1s')),
		libs=m[:,1],
        count=m[:,2],
        fdr=m[:,3],
        pred_fdr=m[:,4],
        odds_ratio=m[:,5]
	)
	return verified_stats
end
function checkinteract(conditions::Vector{Interactions}, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr =0.05)
    verified_stats = DataFrame(name1=String[p[1] for p in verified_pairs], name2=String[p[2] for p in verified_pairs])
    for ints in conditions
        verified_stats = innerjoin(verified_stats, checkinteract(ints, verified_pairs; min_reads=min_reads, max_fdr =max_fdr); on=[:name1, :name2], makeunique=true)
    end
    return verified_stats
end


function uniqueinteract(ints::Interactions; min_reads=5, max_fdr =0.05)
    df = asdataframe(ints; min_reads=min_reads, max_fdr =max_fdr)
    Set(Set((a,b)) for (a,b) in zip(df[!, :name1], df[!, :name2]))
end

function checktops(interact::Interactions; top_cut=20, check=:singles)
    check in (:singles, :interact) || throw(AssertionError("check has to be :singles or :interact"))
    hasfdrvalues(interact) || throw(AssertionError("Assign p-values before running this."))
    top_nodes = checksortperm(interact.nodes[!, check === :singles ? :nb_single : :nb_ints]; rev=true)[1:top_cut]
    filtered_edges = Dict{Tuple{String, String}, Tuple{Float64, Float64}}()
    for edge in eachrow(interact.edges)
        (edge[:src] in top_nodes) && (edge[:dst] in top_nodes) || continue
        push!(filtered_edges, (interact.nodes[edge[:src], :name], interact.nodes[edge[:dst], :name]) => (edge[:fdr], edge[:odds_ratio]))
    end
    return filtered_edges
end

function bp_score_dist_plot(interact::Interactions, genome_model_ecdf::ECDF, randseq_model_ecdf::ECDF, plot_fdr_level::Float64)

    interact_ecdf = ecdf(interact.edges.pred_score[.!isnan.(interact.edges.pred_score)])

    max_score = Int(floor(max(
        maximum(interact_ecdf.sorted_values),
        maximum(randseq_model_ecdf.sorted_values),
        maximum(genome_model_ecdf.sorted_values)
    )))

    nan_index = .!isnan.(interact.edges.pred_fdr)
    si_fdr = sort(collect(enumerate(zip(interact.edges.pred_fdr[nan_index], interact.edges.pred_pvalue[nan_index]))), by=x->x[2][2], rev=true)

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
    sigfish_index = interact.edges.fdr[valid_index] .<= plotting_fdr_level
    h2 = histogram(x=logodds[sigfish_index .& infindex], name="fdr <= $plotting_fdr_level")

    all_index = interact.edges.pred_fdr[valid_index] .<= 1.0
    h3 = histogram(x=logodds[all_index], name="all")
    sigpred_index = interact.edges.pred_fdr[valid_index] .<= plotting_fdr_level
    h4 = histogram(x=logodds[sigpred_index .& infindex], name="fdr <= $plotting_fdr_level")

    plot([h1, h2]), plot([h3, h4])
end

function countdata(data::Vector{Int})
    counts = zeros(Int, maximum(data; init=1))
    for d in data[data .> 0]
        counts[d] += 1
    end
    counts
end

function degreecounts(interact::Interactions; index=trues(nrow(interact.edges)))
    src, dst = interact.edges.src[index], interact.edges.dst[index]
    degs = [length(union!(Set(dst[src .== i]), Set(src[dst .== i]))) for i in 1:nrow(interact.nodes)]
    countdata(degs)
end

function degrees(interact::Interactions; index=trues(nrow(interact.edges)))
    src, dst = interact.edges.src[index], interact.edges.dst[index]
    [length(union!(Set(dst[src .== i]), Set(src[dst .== i]))) for i in 1:nrow(interact.nodes)]
end

function degree_dist_plots(interact::Interactions, plotting_fdr_level::Float64)

    sigfish_index = interact.edges.fdr .<= plotting_fdr_level
    degs = log.(degrees(interact; index=sigfish_index) .+ 1)
    h2 = histogram(x=degs, name="fdr <= $plotting_fdr_level")

    sigpred_index = interact.edges.pred_fdr .<= plotting_fdr_level
    degs = log.(degrees(interact; index=sigpred_index) .+ 1)
    h4 = histogram(x=degs, name="fdr <= $plotting_fdr_level")

    plot(h2), plot(h4)
end

function annotation_type_heatmap(interact::Interactions, plotting_fdr_level::Float64)
    
    type1 = interact.nodes.type[interact.edges.src]
    type2 = interact.nodes.type[interact.edges.dst]
    types = collect(union!(Set(type1), Set(type2)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))

    types_counter = zeros(Float64, (length(types), length(types)))
    type1 = interact.nodes.type[interact.edges.src[interact.edges.fdr .<= plotting_fdr_level]]
    type2 = interact.nodes.type[interact.edges.dst[interact.edges.fdr .<= plotting_fdr_level]]
    for (t1, t2) in zip(type1, type2)
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    types_counter ./= sum(types_counter)
    h1 = heatmap(x=types, y=types, z=copy(types_counter))

    types_counter = zeros(Float64, (length(types), length(types)))
    type1 = interact.nodes.type[interact.edges.src[interact.edges.pred_fdr .<= plotting_fdr_level]]
    type2 = interact.nodes.type[interact.edges.dst[interact.edges.pred_fdr .<= plotting_fdr_level]]
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
function basepairing_string(aln::PairwiseAlignment, offset1::Int, offset2::Int)
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    posw = ndigits(max(abs(offset1 + anchors[end].seqpos), abs(offset2 - anchors[1].refpos))) + 2

    i = 0
    seqpos = offset1 + anchors[1].seqpos
    refpos = offset2 - anchors[1].refpos + 2
    seqbuf = IOBuffer()
    refbuf = IOBuffer()
    matbuf = IOBuffer()

    print(seqbuf, "RNA1:", lpad(seqpos>0 ? "+$seqpos" : "$(seqpos-1)", posw), ' ')
    print(refbuf, "RNA2:", lpad(refpos>0 ? "+$refpos" : "$(refpos-1)", posw), ' ')
    print(matbuf, " "^(posw + 6))

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
    print(matbuf)

    String(take!(seqbuf)) * "<br>" * String(take!(matbuf)) * "<br>" * String(take!(refbuf))
end
function alignment_ascii_plot(i1::Int, i2::Int, p1::Int, p2::Int, interact::Interactions, genome::Dict{String, BioSequences.LongDNA{4}}, check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel)

    ref1::String, strand1::Char, l1::Int, r1::Int, c1::Int = interact.nodes[i1, [:ref, :strand, :left, :right, :cds]]
    ref2::String, strand2::Char, l2::Int, r2::Int, c2::Int = interact.nodes[i2, [:ref, :strand, :left, :right, :cds]]

    s1 = strand1=='+' ?
        genome[ref1][(p1-check_interaction_distances[1]):(p1-check_interaction_distances[2])] :
        BioSequences.reverse_complement(genome[ref1][(p1+check_interaction_distances[2]):(p1+check_interaction_distances[1])])
    s2 = strand2=='-' ?
        BioSequences.complement(genome[ref2][(p2-check_interaction_distances[1]):(p2-check_interaction_distances[2])]) :
        BioSequences.reverse(genome[ref2][(p2+check_interaction_distances[2]):(p2+check_interaction_distances[1])])
    p = pairalign(LocalAlignment(), s1, s2, model)
    basepairing_string(alignment(p),
        (strand1=='+' ? ((p1-check_interaction_distances[1])-(c1>0 ? c1 : l1)) : ((c1>0 ? c1 : r1)-(p1+check_interaction_distances[1]))),
        (strand2=='-' ? ((c2>0 ? c2 : r2)-(p2-check_interaction_distances[1])) : ((p2+check_interaction_distances[1])-(c2>0 ? c2 : l2))))
end