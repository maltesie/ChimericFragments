hasfdrvalues(ints::Interactions) = "fdr" in names(ints.edges)

function checkinteractions(ints::Interactions, verified_pairs::Vector{Tuple{String,String}};
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
function checkinteractions(conditions::Vector{Interactions}, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr =0.05)
    verified_stats = DataFrame(name1=String[p[1] for p in verified_pairs], name2=String[p[2] for p in verified_pairs])
    for ints in conditions
        verified_stats = innerjoin(verified_stats, checkinteractions(ints, verified_pairs; min_reads=min_reads, max_fdr =max_fdr); on=[:name1, :name2], makeunique=true)
    end
    return verified_stats
end
function checkinteractions(interaction_files::SingleTypeFiles, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr =0.05)
    interaction_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(interaction_file) for interaction_file in interaction_files]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr =max_fdr)
end
function checkinteractions(interaction_files::SingleTypeFiles, verified_pairs_file::String; min_reads=5, max_fdr =0.05)
    interaction_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(interaction_file) for interaction_file in interaction_files]
    df = DataFrame(CSV.File(verified_pairs_file; stringtype=String))
    verified_pairs = [(a,b) for (a,b) in eachrow(df)]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr =max_fdr)
end

function uniqueinteractions(ints::Interactions; min_reads=5, max_fdr =0.05)
    df = asdataframe(ints; min_reads=min_reads, max_fdr =max_fdr)
    Set(Set((a,b)) for (a,b) in zip(df[!, :name1], df[!, :name2]))
end

function checktops(interactions::Interactions; top_cut=20, check=:singles)
    check in (:singles, :interactions) || throw(AssertionError("check has to be :singles or :interactions"))
    hasfdrvalues(interactions) || throw(AssertionError("Assign p-values before running this."))
    top_nodes = checksortperm(interactions.nodes[!, check === :singles ? :nb_single : :nb_ints]; rev=true)[1:top_cut]
    filtered_edges = Dict{Tuple{String, String}, Tuple{Float64, Float64}}()
    for edge in eachrow(interactions.edges)
        (edge[:src] in top_nodes) && (edge[:dst] in top_nodes) || continue
        push!(filtered_edges, (interactions.nodes[edge[:src], :name], interactions.nodes[edge[:dst], :name]) => (edge[:fdr], edge[:odds_ratio]))
    end
    return filtered_edges
end

function bp_score_dist_plot(interactions::Interactions, genome_model_ecdf::ECDF, randseq_model_ecdf::ECDF, max_bp_fdr::Float64)

    interactions_ecdf = ecdf(Int.(interactions.edges.pred_score[.!isnan.(interactions.edges.pred_score)]))

    max_score = Int(max(
        maximum(interactions_ecdf.sorted_values),
        maximum(randseq_model_ecdf.sorted_values),
        maximum(genome_model_ecdf.sorted_values)
    ))

    nan_index = .!isnan.(interactions.edges.pred_fdr)
    si_fdr = sort(collect(enumerate(zip(interactions.edges.pred_fdr[nan_index], interactions.edges.pred_pvalue[nan_index]))), by=x->x[2][2], rev=true)

    fdr01_p_index = findfirst(x->x[2][1]<=max_bp_fdr, si_fdr)
    fdr01_p = isnothing(fdr01_p_index) ? 1.0 : si_fdr[fdr01_p_index][2][2]
    genome_model_fdr_01_score = findfirst(x->((1-genome_model_ecdf(x))<fdr01_p), 1:(max_score+1)) - 1

    genome_model_pdf = diff(genome_model_ecdf.(1:(max_score+1)))
    randseq_model_pdf = diff(randseq_model_ecdf.(1:(max_score+1)))
    interactions_pdf = diff(interactions_ecdf.(1:(max_score+1)))

    p = plot(1:max_score, randseq_model_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 1, label = "random", legend = :topright)
    plot!(p, 1:max_score, genome_model_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 2, label = "random, from genome")
    plot!(p, 1:max_score, interactions_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 3, label = "around ligation points")
    vline!(p, [genome_model_fdr_01_score], label="fdr = $max_bp_fdr")
    return p
end

function alignment_histogram(l1::Vector{Float64}, r1::Vector{Float64}, l2::Vector{Float64}, r2::Vector{Float64},
            bins::AbstractRange, max_fdr::Float64, fdrs::Vector{Float64})

    sig_index = fdrs .<= max_fdr
    nonsig_index = .!sig_index

    h1 = histogram(l1; bins=bins, label="all", legend=:topright)
    histogram!(h1, l1[nonsig_index]; bins=bins, label="fdr > $max_fdr")
    histogram!(h1, l1[sig_index]; bins=bins, label="fdr <= $max_fdr")
    title!(h1, "RNA1, left end")

    h2 = histogram(r1; bins=bins, label="all", legend=:topleft)
    histogram!(h2, r1[nonsig_index]; bins=bins, label="fdr > $max_fdr")
    histogram!(h2, r1[sig_index]; bins=bins, label="fdr <= $max_fdr")
    title!(h2, "RNA1, right end")

    h3 = histogram(l2; bins=bins, label="all", legend=:topright)
    histogram!(h3, l2[nonsig_index]; bins=bins, label="fdr > $max_fdr")
    histogram!(h3, l2[sig_index]; bins=bins, label="fdr <= $max_fdr")
    title!(h3, "RNA2, left end")

    h4 = histogram(r2; bins=bins, label="all", legend=:topleft)
    histogram!(h4, r2[nonsig_index]; bins=bins, label="fdr > $max_fdr")
    histogram!(h4, r2[sig_index]; bins=bins, label="fdr <= $max_fdr")
    title!(h4, "RNA2, right end")

    h5 = histogram(r1 .- l1 .+ 1; bins=bins, label="all", legend=:topright)
    histogram!(h5, r1[nonsig_index] .- l1[nonsig_index] .+ 1; bins=bins, label="fdr > $max_fdr")
    histogram!(h5, r1[sig_index] .- l1[sig_index] .+ 1; bins=bins, label="fdr <= $max_fdr")
    title!(h5, "RNA1, length")

    h6 = histogram(r2 .- l2 .+ 1; bins=bins, label="all", legend=:topright)
    histogram!(h6, r2[nonsig_index] .- l2[nonsig_index] .+ 1; bins=bins, label="fdr > $max_fdr")
    histogram!(h6, r2[sig_index] .- l2[sig_index] .+ 1; bins=bins, label="fdr <= $max_fdr")
    title!(h6, "RNA2, length")

    plot(h1, h5, h2, h3, h6, h4; layout=(2,3), size=(1800,800))
end

function bp_clipping_dist_plots(interactions::Interactions, bp_distance::Tuple{Int,Int}, max_fisher_fdr::Float64, max_bp_fdr::Float64)

    bins = 0:(bp_distance[1]-bp_distance[2])
    nan_index = .!isnan.(interactions.edges.pred_fdr)

    l1 = interactions.edges.pred_cl1[nan_index]
    r1 = interactions.edges.pred_cr1[nan_index]
    l2 = interactions.edges.pred_cl2[nan_index]
    r2 = interactions.edges.pred_cr2[nan_index]

    bp_fdrs = interactions.edges.pred_fdr[nan_index]
    p1 = alignment_histogram(l1, r1, l2, r2, bins, max_bp_fdr, bp_fdrs)

    fisher_fdrs = interactions.edges.fdr[nan_index]
    p2 = alignment_histogram(l1, r1, l2, r2, bins, max_fisher_fdr, fisher_fdrs)

    return p1, p2
end

function log_chimeric_counts_distribution_plot(interactions::Interactions, max_fisher_fdr::Float64, max_bp_fdr::Float64)
    counts = log10.(interactions.edges.nb_ints)
    h = histogram(counts; label="all", legend=:topright)
    title!("counts")
    sigfish_index = interactions.edges.fdr .<= max_fisher_fdr
    histogram!(h, counts[sigfish_index]; label="fisher fdr <= $max_fisher_fdr", legend=:topright)
    sigpred_index = interactions.edges.pred_fdr .<= max_bp_fdr
    histogram!(h, counts[sigpred_index]; label="bp fdr <= $max_bp_fdr", legend=:topright)
end

function log_odds_ratio_distribution_plot(interactions::Interactions, max_fisher_fdr::Float64, max_bp_fdr::Float64)
    odds = log.(interactions.edges.odds_ratio)
    h = histogram(odds; label="all", legend=:topright)
    title!("log odds ratio")
    sigfish_index = interactions.edges.fdr .<= max_fisher_fdr
    histogram!(h, odds[sigfish_index]; label="fisher fdr <= $max_fisher_fdr", legend=:topright)
    sigpred_index = interactions.edges.pred_fdr .<= max_bp_fdr
    histogram!(h, odds[sigpred_index]; label="bp fdr <= $max_bp_fdr", legend=:topright)
end

