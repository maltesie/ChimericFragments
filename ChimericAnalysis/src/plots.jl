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

function bp_score_dist_plot(interactions::Interactions, genome_model_ecdf::ECDF, randseq_model_ecdf::ECDF, plot_fdr_levels::Vector{Float64})

    interactions_ecdf = ecdf(Int.(interactions.edges.pred_score[.!isnan.(interactions.edges.pred_score)]))

    max_score = Int(max(
        maximum(interactions_ecdf.sorted_values),
        maximum(randseq_model_ecdf.sorted_values),
        maximum(genome_model_ecdf.sorted_values)
    ))

    nan_index = .!isnan.(interactions.edges.pred_fdr)
    si_fdr = sort(collect(enumerate(zip(interactions.edges.pred_fdr[nan_index], interactions.edges.pred_pvalue[nan_index]))), by=x->x[2][2], rev=true)

    genome_model_pdf = diff(genome_model_ecdf.(1:(max_score+1)))
    randseq_model_pdf = diff(randseq_model_ecdf.(1:(max_score+1)))
    interactions_pdf = diff(interactions_ecdf.(1:(max_score+1)))

    p = plot(1:max_score, randseq_model_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 1, label = "random",
        legend = :topright, size=(600,400))
    xlabel!(p, "affine gap model score")
    ylabel!(p, "empirical density")
    plot!(p, 1:max_score, genome_model_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 2, label = "random, from genome")
    plot!(p, 1:max_score, interactions_pdf, fillrange = zeros(max_score), fillalpha = 0.35, c = 3, label = "around ligation points")

    for max_bp_fdr in plot_fdr_levels
        fdr_p_index = findfirst(x->x[2][1]<=max_bp_fdr, si_fdr)
        fdr_p = isnothing(fdr_p_index) ? 1.0 : si_fdr[fdr_p_index][2][2]
        genome_model_fdr_score = findfirst(x->((1-genome_model_ecdf(x))<fdr_p), 1:(max_score+1)) - 1
        vline!(p, [genome_model_fdr_score], label="fdr = $max_bp_fdr")
    end

    return p
end

function alignment_histogram(l1::Vector{Float64}, r1::Vector{Float64}, l2::Vector{Float64}, r2::Vector{Float64},
            bins::AbstractRange, max_fdrs::Vector{Float64}, fdrs::Vector{Float64})

    h1 = histogram(l1; bins=bins, label="all", legend=:topright)
    title!(h1, "RNA1, left end")
    ylabel!(h1, "count")
    h2 = histogram(r1; bins=bins, label="all", legend=:topleft)
    title!(h2, "RNA1, right end")
    h3 = histogram(l2; bins=bins, label="all", legend=:topright)
    ylabel!(h3, "count")
    xlabel!(h3, "position in alignment")
    title!(h3, "RNA2, left end")
    h4 = histogram(r2; bins=bins, label="all", legend=:topleft)
    xlabel!(h4, "position in alignment")
    title!(h4, "RNA2, right end")
    h5 = histogram(r1 .- l1 .+ 1; bins=bins, label="all", legend=:topright)
    title!(h5, "RNA1, length")
    h6 = histogram(r2 .- l2 .+ 1; bins=bins, label="all", legend=:topright)
    xlabel!(h6, "position in alignment")
    title!(h6, "RNA2, length")

    for max_fdr in reverse(max_fdrs)
        sig_index = fdrs .<= max_fdr
        histogram!(h1, l1[sig_index]; bins=bins, label="fdr <= $max_fdr")
        histogram!(h2, r1[sig_index]; bins=bins, label="fdr <= $max_fdr")
        histogram!(h3, l2[sig_index]; bins=bins, label="fdr <= $max_fdr")
        histogram!(h4, r2[sig_index]; bins=bins, label="fdr <= $max_fdr")
        histogram!(h5, r1[sig_index] .- l1[sig_index] .+ 1; bins=bins, label="fdr <= $max_fdr")
        histogram!(h6, r2[sig_index] .- l2[sig_index] .+ 1; bins=bins, label="fdr <= $max_fdr")
    end

    plot(h1, h5, h2, h3, h6, h4; layout=(2,3), size=(1800,800), margin=7mm)
end

function bp_clipping_dist_plots(interactions::Interactions, bp_distance::Tuple{Int,Int}, plotting_fdr_levels::Vector{Float64})

    bins = 0:(bp_distance[1]-bp_distance[2])
    nan_index = .!isnan.(interactions.edges.pred_fdr)

    l1 = interactions.edges.pred_cl1[nan_index]
    r1 = interactions.edges.pred_cr1[nan_index]
    l2 = interactions.edges.pred_cl2[nan_index]
    r2 = interactions.edges.pred_cr2[nan_index]

    bp_fdrs = interactions.edges.pred_fdr[nan_index]
    p1 = alignment_histogram(l1, r1, l2, r2, bins, plotting_fdr_levels, bp_fdrs)

    fisher_fdrs = interactions.edges.fdr[nan_index]
    p2 = alignment_histogram(l1, r1, l2, r2, bins, plotting_fdr_levels, fisher_fdrs)

    return p1, p2
end

function countdata(data::Vector{Int})
    counts = zeros(Int, maximum(data; init=1))
    for d in data[data .> 0]
        counts[d] += 1
    end
    counts
end

function fisher_pred_histograms(data::Vector{Float64}, fisher_fdrs::Vector{Float64}, bp_fdrs::Vector{Float64},
        bins::AbstractRange, plotting_fdr_levels::Vector{Float64}, title::String)

    h1 = histogram(data; label="all", bins=bins, legend=:topright)
    title!(h1, "fisher fdr")
    for max_fisher_fdr in reverse(plotting_fdr_levels)
        sigfish_index = fisher_fdrs .<= max_fisher_fdr
        histogram!(h1, data[sigfish_index]; label="fisher fdr <= $max_fisher_fdr", bins=bins)
    end

    nan_index = .!isnan.(bp_fdrs)
    h2 = histogram(data[nan_index]; label="all", bins=bins, legend=:topright)
    title!(h2, "basepairing prediction fdr")
    for max_bp_fdr in reverse(plotting_fdr_levels)
        sigpred_index = bp_fdrs .<= max_bp_fdr
        histogram!(h2, data[sigpred_index]; label="bp fdr <= $max_bp_fdr", bins=bins)
    end

    plot(h1, h2; layout=(1,2), size=(1200,400), margin=5mm, xlabel=title, ylabel="count")
end

function fisher_pred_lineplots(data::Vector{Int}, fisher_fdrs::Vector{Float64}, bp_fdrs::Vector{Float64}, plotting_fdr_levels::Vector{Float64}, title::String)

    counts = countdata(data)
    p1 = plot(counts; label="all", xaxis=:log, fillalpha=0.35, legend=:topright)
    title!(p1, "fisher fdr")
    for max_fisher_fdr in reverse(plotting_fdr_levels)
        sigfish_index = fisher_fdrs .<= max_fisher_fdr
        counts = countdata(data[sigfish_index])
        plot!(p1, counts; label="fisher fdr <= $max_fisher_fdr", xaxis=:log, fillalpha=0.35)
    end

    nan_index = .!isnan.(bp_fdrs)
    counts = countdata(data[nan_index])
    p2 = plot(counts; label="all", xaxis=:log, fillalpha=0.35, legend=:topright)
    title!(p2, "basepairing prediction fdr")
    for max_bp_fdr in reverse(plotting_fdr_levels)
        sigpred_index = bp_fdrs .<= max_bp_fdr
        counts = countdata(data[sigpred_index])
        plot!(p2, counts; label="bp fdr <= $max_bp_fdr", xaxis=:log, fillalpha=0.35)
    end

    plot(p1, p2; layout=(1,2), size=(1200,400), margin=5mm, xlabel=title, ylabel="count")
end

function interaction_distribution_plots(interactions::Interactions, plotting_fdr_levels::Vector{Float64})

    counts = interactions.edges.nb_ints
    p1 = fisher_pred_lineplots(counts, interactions.edges.fdr, interactions.edges.pred_fdr, plotting_fdr_levels, "log read count")

    odds = log.(interactions.edges.odds_ratio)
    inf_index = .!isinf.(odds)
    bins = sum(inf_index) > 0 ? range(minimum(odds[inf_index]), maximum(odds[inf_index]), 100) : range(-1,1,100)
    p2 = fisher_pred_histograms(odds, interactions.edges.fdr, interactions.edges.pred_fdr, bins, plotting_fdr_levels, "log odds ratio")

    return p1, p2
end

function degreecounts(interactions::Interactions; index=trues(nrow(interactions.edges)))
    src, dst = interactions.edges.src[index], interactions.edges.dst[index]
    degs = [length(union!(Set(dst[src .== i]), Set(dst[src .== i]))) for i in 1:nrow(interactions.nodes)]
    countdata(degs)
end

function node_distribution_plot(interactions::Interactions, plotting_fdr_levels::Vector{Float64})

    counts = degreecounts(interactions)
    p1 = plot(counts; label="all", xaxis=:log, fillalpha=0.35, legend=:topright)
    title!(p1, "fisher fdr")

    for max_fisher_fdr in reverse(plotting_fdr_levels)
        sigfish_index = interactions.edges.fdr .<= max_fisher_fdr
        counts = degreecounts(interactions; index=sigfish_index)
        plot!(p1, counts; label="fisher fdr <= $max_fisher_fdr", xaxis=:log, fillalpha=0.35)
    end

    nan_index = .!isnan.(interactions.edges.pred_fdr)
    counts = degreecounts(interactions; index=nan_index)
    p2 = plot(counts; label="all", xaxis=:log, fillalpha=0.35, legend=:topright)
    title!(p2, "basepairing prediction fdr")

    for max_bp_fdr in reverse(plotting_fdr_levels)
        sigpred_index = interactions.edges.pred_fdr .<= max_bp_fdr
        counts = degreecounts(interactions; index=sigpred_index)
        plot!(p2, counts; label="bp fdr <= $max_bp_fdr", xaxis=:log, fillalpha=0.35)
    end

    plot(p1, p2; layout=(1,2), size=(1200,400), margin=6mm, xlabel="log degree", ylabel="count")
end

function annotation_type_heatmap(interactions::Interactions, plotting_fdr_levels::Vector{Float64})
    type1 = interactions.nodes.type[interactions.edges.src]
    type2 = interactions.nodes.type[interactions.edges.dst]
    types = collect(union!(Set(type1), Set(type2)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))

    types_counter = zeros(Int, (length(types), length(types)))
    for (t1, t2) in zip(type1, type2)
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    plots = [heatmap(types, types, types_counter; legend=:none)]
    titles = ["all interactions"]
    for pl in plotting_fdr_levels
        types_counter = zeros(Int, (length(types), length(types)))
        type1 = interactions.nodes.type[interactions.edges.src[interactions.edges.fdr .<= pl]]
        type2 = interactions.nodes.type[interactions.edges.dst[interactions.edges.fdr .<= pl]]
        for (t1, t2) in zip(type1, type2)
            types_counter[type_trans[t1], type_trans[t2]] += 1
        end
        push!(plots, heatmap(types, types, types_counter; legend=:none))
        push!(titles, "max fisher fdr = $pl")
    end
    for pl in plotting_fdr_levels
        types_counter = zeros(Int, (length(types), length(types)))
        type1 = interactions.nodes.type[interactions.edges.src[interactions.edges.pred_fdr .<= pl]]
        type2 = interactions.nodes.type[interactions.edges.dst[interactions.edges.pred_fdr .<= pl]]
        for (t1, t2) in zip(type1, type2)
            types_counter[type_trans[t1], type_trans[t2]] += 1
        end
        push!(plots, heatmap(types, types, types_counter; legend=:none))
        push!(titles, "max bp fdr = $pl")
    end
    l = @layout[a{0.25w} grid(2,length(plotting_fdr_levels))]
    plot(plots...; layout=l, title=reshape(titles, (1, length(titles))), size=(600*(length(plotting_fdr_levels)+1),800), margin=9mm)
end