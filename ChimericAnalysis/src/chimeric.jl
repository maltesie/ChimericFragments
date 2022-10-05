struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_ints=>Int[], :nb_ints_src=>Int[], :nb_ints_dst=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :meanlength1=>Float64[], :meanlength2=>Float64[], :nms1=>Float64[], :nms2=>Float64[])
    edgestats = Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}()
    Interactions(nodes, edges, edgestats, Symbol[])
end

Interactions(alignments::Alignments; replicate_id=:first, min_distance=1000, max_ligation_distance=5, filter_types=[], multi_detection_method=:annotation) =
    append!(Interactions(), alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                filter_types=filter_types, multi_detection_method=multi_detection_method)

"""
Load Interactions struct from jld2 file.
"""
Interactions(filepath::String) = load(filepath, "interactions")

Base.length(interactions::Interactions) = nrow(interactions.edges)

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

struct MergedAlignedPart
    ref::Interval{AlignmentAnnotation}
    seq1::UnitRange{Int}
    seq2::UnitRange{Int}
    nms1::UInt32
    nms2::UInt32
    read::Symbol
end

MergedAlignedPart(part::AlignedPart) = MergedAlignedPart(part.ref, isfirstread(part) ? part.seq : (1:0), isfirstread(part) ? (1:0) : part.seq, part.nms, part.nms, part.read)
function MergedAlignedPart(part1::AlignedPart, part2::AlignedPart)
    firstpart, secondpart = isfirstread(part1) ? (part1, part2) : (part2, part1)
    ref = Interval(refname(firstpart), min(leftposition(firstpart), leftposition(secondpart)),
        max(rightposition(firstpart), rightposition(secondpart)), strand(firstpart),
        annotation(overlap(firstpart)>overlap(secondpart) ? firstpart : secondpart))
    return MergedAlignedPart(ref, firstpart.seq, secondpart.seq, firstpart.nms, secondpart.nms, :both)
end

function mergedparts(parts::Vector{AlignedPart})
    merged_parts = MergedAlignedPart[]
    already_merged = zeros(Bool, length(parts))
    for (i, part) in enumerate(parts)
        already_merged[i] && continue
        nextolpi = findnext(x->isoverlapping(parts[i].ref, x.ref) && !sameread(parts[i], x), parts, i+1)
        if isnothing(nextolpi)
            push!(merged_parts, MergedAlignedPart(part))
        else
            push!(merged_parts, MergedAlignedPart(part, parts[nextolpi]))
            already_merged[nextolpi] = true
        end
    end
    return merged_parts
end

isligationpoint(part1::MergedAlignedPart, part2::MergedAlignedPart; max_distance=5) =
    ((isempty(part1.seq1) || isempty(part2.seq1)) ? false : (first(part2.seq1) - last(part1.seq1) <= max_distance)) ||
    ((isempty(part1.seq2) || isempty(part2.seq2)) ? false : (first(part2.seq2) - last(part1.seq2) <= max_distance))

Base.length(p::MergedAlignedPart) = rightposition(p.ref)-leftposition(p.ref)+1

function ischimeric(part1::MergedAlignedPart, part2::MergedAlignedPart; min_distance=1000, check_annotation=true)
    check_annotation && (part1.ref.metadata.name == part2.ref.metadata.name) && return false
    return distance(part1.ref, part2.ref) > min_distance
end

function myhash(part::AlignedPart; use_type=true)
    return use_type ? hash(name(part)*type(part)) : hash(name(part))
end

function myhash(part::MergedAlignedPart; use_type=true)
    return use_type ? hash(part.ref.metadata.name*part.ref.metadata.type) : hash(part.ref.metadata.name)
end

function Base.append!(interactions::Interactions, alignments::Alignments, replicate_id::Symbol;
                        min_distance=1000, max_ligation_distance=5, filter_types=[], multi_detection_method=:annotation)
    if !(String(replicate_id) in interactions.replicate_ids)
        interactions.edges[:, replicate_id] = zeros(Int, nrow(interactions.edges))
        push!(interactions.replicate_ids, replicate_id)
    end
    trans = Dict{UInt, Int}(interactions.nodes[i, :hash]=>i for i in 1:nrow(interactions.nodes))
    edgestats = interactions.edgestats
    exclude_count = 0
    single_count = 0
    chimeric_count = 0
    total_count = 0
    for alignment in alignments
        total_count += 1
        !isempty(filter_types) && typein(alignment, filter_types) && (exclude_count += 1; continue)
        is_chimeric = ischimeric(alignment; min_distance=min_distance)
        is_chimeric ? chimeric_count += 1 : single_count +=1
        is_multi = is_chimeric ? ismulti(alignment; method=multi_detection_method) : false
        alnparts = parts(alignment)
        filter!(x->hasannotation(x), alnparts)
        merged_parts = mergedparts(alnparts)

        for part in alnparts
            h = myhash(part)
            h in keys(trans) && continue
            if !(h in keys(trans))
                trans[h] = length(trans) + 1
                push!(interactions.nodes, (name(part), type(part), refname(part), 0, 0, 0, 0, strand(part), h))
            end
            is_chimeric || (interactions.nodes[trans[h], :nb_single] += 1)
        end

        for (part1, part2) in combinations(merged_parts, 2)
            ischimeric(part1, part2; min_distance=min_distance) || continue
            a, b = trans[myhash(part1)], trans[myhash(part2)]
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
            leftpos = rightposition(part1.ref)
            rightpos = leftposition(part2.ref)
            leftcounter, rightcounter = isligationpoint(part1, part2; max_distance=max_ligation_distance) ? (leftligationcounter, rightligationcounter) : (leftintcounter, rightintcounter)
            leftpos in keys(leftcounter) ? (leftcounter[leftpos]+=1) : (leftcounter[leftpos]=1)
            rightpos in keys(rightcounter) ? (rightcounter[rightpos]+=1) : (rightcounter[rightpos]=1)
            for (s,v) in zip((:meanlength1, :meanlength2, :nms1, :nms2),(length(part1), length(part2), max(part1.nms1, part1.nms2), max(part2.nms1, part2.nms2)))
                interactions.edges[iindex, s] += (v - interactions.edges[iindex, s]) / interactions.edges[iindex, :nb_ints]
            end
        end
    end
    println("Processed $total_count reads, found $chimeric_count chimeras and $single_count singles and excluded $exclude_count.")
    return interactions
end

hasfdrvalues(ints::Interactions) = "fdr" in names(ints.edges)

function checkinteractions(ints::Interactions, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05, check_uppercase=true)
	verified_dict = merge(Dict(pair=>zeros(2) for pair in verified_pairs),
							Dict(reverse(pair)=>zeros(2) for pair in verified_pairs))
    check_uppercase && (verified_dict = Dict((uppercase(p[1]), uppercase(p[2]))=>v for (p,v) in verified_dict))
	df = asdataframe(ints; min_reads=min_reads, max_fdr=max_fdr)
    for row in eachrow(df)
        key = (row[:name1], row[:name2])
        check_uppercase && (key = (uppercase(key[1]), uppercase(key[2])))
        if key in keys(verified_dict)
            verified_dict[key][1] = row[:in_libs]
            verified_dict[key][2] = row[:nb_ints]
        end
    end

	sorted_keys = vcat([[pair, reverse(pair)] for pair in verified_pairs]...)
    check_uppercase && (sorted_keys = [(uppercase(p[1]), uppercase(p[2])) for p in sorted_keys])
	m = reduce(hcat, [verified_dict[key] for key in sorted_keys])'
	verified_stats = DataFrame(
		name1=String[n[1] for n in verified_pairs],
		name2=String[n[2] for n in verified_pairs],
		libs=max.(m[:,1][1:2:end], m[:,1][2:2:end]), count=m[:,2][1:2:end] .+ m[:,2][2:2:end]
		)
	return verified_stats
end
function checkinteractions(conditions::Vector{Interactions}, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05)
    verified_stats = DataFrame(name1=String[p[1] for p in verified_pairs], name2=String[p[2] for p in verified_pairs])
    for ints in conditions
        verified_stats = innerjoin(verified_stats, checkinteractions(ints, verified_pairs; min_reads=min_reads, max_fdr=max_fdr); on=[:name1, :name2], makeunique=true)
    end
    return verified_stats
end
function checkinteractions(interaction_files::SingleTypeFiles, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05)
    interaction_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(interaction_file) for interaction_file in interaction_files]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr=max_fdr)
end
function checkinteractions(interaction_files::SingleTypeFiles, verified_pairs_file::String; min_reads=5, max_fdr=0.05)
    interaction_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(interaction_file) for interaction_file in interaction_files]
    df = DataFrame(CSV.File(verified_pairs_file; stringtype=String))
    verified_pairs = [(a,b) for (a,b) in eachrow(df)]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr=max_fdr)
end

function uniqueinteractions(ints::Interactions; min_reads=5, max_fdr=0.05)
    df = asdataframe(ints; min_reads=min_reads, max_fdr=max_fdr)
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

function addpvalues!(interactions::Interactions; method=:fisher, fisher_tail=:right, include_read_identity=true, include_singles=true)
    @assert method in (:disparity, :fisher)
    pvalues = ones(nrow(interactions.edges))

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

    include_singles && (other_source .+= interactions.nodes[interactions.edges[!, :src], :nb_single])
    include_singles && (other_target .+= interactions.nodes[interactions.edges[!, :dst], :nb_single])
    total_other = sum(interactions.edges[!, :nb_ints]) .- ints_between .- other_source .- other_target .+ (include_singles ? sum(interactions.nodes[!, :nb_single]) : 0)

    odds_ratio = (ints_between .* total_other) ./ (other_target .* other_source)

    if method === :fisher
        tests = FisherExactTest.(ints_between, other_target, other_source, total_other)
        pvalues = pvalue.(tests; tail=fisher_tail)
    else
        throw(AssertionError("$method not supported!"))
    end

    adjp = adjust(PValues(pvalues), BenjaminiHochberg())
    interactions.edges[:, :odds_ratio] = odds_ratio
    interactions.edges[:, :p_value] = pvalues
    interactions.edges[:, :fdr] = adjp
    return interactions
end

function addpositions!(interactions::Interactions, features::Features)
    tus = Dict(hash(name(feature)*type(feature))=>(leftposition(feature), rightposition(feature)) for feature in features)
    for (col_int, col_float) in zip([:left1, :right1, :left2, :right2], [:relposints1, :relposints2, :relposligation1, :relposligation2])
        interactions.edges[:, col_int] = Vector{Int}(undef, length(interactions))
        interactions.edges[:, col_float] = Vector{Float64}(undef, length(interactions))
    end
    for edge_row in eachrow(interactions.edges)
        (feature1_left, feature1_right) = tus[interactions.nodes[edge_row[:src], :hash]]
        (feature2_left, feature2_right) = tus[interactions.nodes[edge_row[:dst], :hash]]
        isnegative1 = interactions.nodes[edge_row[:src], :strand] === '-'
        isnegative2 = interactions.nodes[edge_row[:dst], :strand] === '-'
        stats = interactions.edgestats[(edge_row[:src], edge_row[:dst])]
        p1 = (length(stats[2]) > 0 ? sum(v*p for (v,p) in stats[2]) : Inf)/edge_row.nb_ints, (length(stats[3]) > 0 ? sum(v*p for (v,p) in stats[3]) : Inf)/edge_row.nb_ints
        p2 = (length(stats[4]) > 0 ? sum(v*p for (v,p) in stats[2]) : Inf)/edge_row.nb_ints, (length(stats[3]) > 0 ? sum(v*p for (v,p) in stats[5]) : Inf)/edge_row.nb_ints
        (relposints1, relposligation1) = (p1 .- feature1_left) ./ (feature1_right - feature1_left)
        (relposints2, relposligation2) = (p2 .- feature2_left) ./ (feature2_right - feature2_left)
        isnegative1 && ((relposints1, relposligation1) = (1-relposints1, 1-relposligation1))
        isnegative2 && ((relposints2, relposligation2) = (1-relposints2, 1-relposligation2))
        left1 = minimum(keys(merge(stats[2],stats[3])))
        right1 = maximum(keys(merge(stats[2],stats[3])))
        left2 = minimum(keys(merge(stats[4],stats[5])))
        right2 = maximum(keys(merge(stats[4],stats[5])))
        edge_row[[:left1, :right1, :left2, :right2, :relposints1, :relposints2, :relposligation1, :relposligation2]] =
            round.((left1, right1, left2, right2, relposints1, relposints2, relposligation1, relposligation2); digits=4)
    end
    return interactions
end

function asdataframe(interactions::Interactions; output=:edges, min_reads=5, max_fdr=0.05, include_pvalues=true, include_relative_positions=true, hist_bins=100)
    out_df = copy(interactions.edges)
    filter_index = (out_df[!, :nb_ints] .>= min_reads)
    "fdr" in names(out_df) && (filter_index .&= (out_df[!, :fdr] .<= max_fdr))
    out_df = out_df[filter_index, :]
    if output === :edges
        out_df[!, :meanlength1] = Int.(round.(out_df[!, :meanlength1]))
        out_df[!, :meanlength2] = Int.(round.(out_df[!, :meanlength2]))
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
        out_columns = [:name1, :type1, :ref1, :name2, :type2, :ref2, :nb_ints, :nb_multi, :in_libs, :strand1,
        :left1, :right1, :left2, :right2, :strand2, :meanlength1, :meanlength2, :nms1, :nms2]
        include_pvalues && (out_columns = [out_columns[1:9]..., :p_value, :fdr, out_columns[10:end]...])
        include_relative_positions && (out_columns = [out_columns..., :relposints1, :relposints2, :relposligation1, :relposligation2])
        return sort(out_df[!, out_columns], :nb_ints; rev=true)
    elseif output === :nodes
        out_nodes = copy(interactions.nodes)
        for (i,row) in enumerate(eachrow(out_nodes))
            row[:nb_ints] = sum(out_df[(out_df.src .== i) .| (out_df.dst .== i), :nb_ints])
        end
        return sort(out_nodes[!, [:name, :type, :ref, :nb_single, :nb_ints]], :nb_single; rev=true)
    elseif output === :stats
        stats_df = DataFrame(zeros(nrow(out_df), 4*hist_bins), [
            ["$(i)_ints1" for i in 1:hist_bins]...,
            ["$(i)_ints2" for i in 1:hist_bins]...,
            ["$(i)_ligation1" for i in 1:hist_bins]...,
            ["$(i)_ligation2" for i in 1:hist_bins]...
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
        return sort(stats_df, :nb_ints; rev=true)
    else
        throw(AssertionError("output has to be one of :edges, :nodes, :stats"))
    end
end

function chimeric_analysis(features::Features, bams::SingleTypeFiles, results_path::String, conditions::Dict{String, Vector{Int}};
                            filter_types=["rRNA", "tRNA"], min_distance=1000, prioritize_type="sRNA", min_prioritize_overlap=0.8,
                            overwrite_type="IGR", max_ligation_distance=5, is_reverse_complement=true,
                            include_secondary_alignments=true, include_alternative_alignments=false, model=:fisher, min_reads=5, max_fdr=0.05,
                            overwrite_existing=false, include_read_identity=true, include_singles=true, multi_detection_method=:annotation)

    isdir(joinpath(results_path, "interactions")) || mkpath(joinpath(results_path, "interactions"))
    isdir(joinpath(results_path, "stats")) || mkpath(joinpath(results_path, "stats"))
    isdir(joinpath(results_path, "singles")) || mkpath(joinpath(results_path, "singles"))
    for (condition, r) in conditions
        println("Collecting $(length(r)) samples for condition $condition:")
        if !overwrite_existing &&
            isfile(joinpath(results_path, "interactions", "$(condition).csv")) &&
            isfile(joinpath(results_path, "singles", "$(condition).csv")) &&
            isfile(joinpath(results_path, "stats", "$(condition).csv"))
            println("Found results files. Skipping...")
            continue
        end
        replicate_ids = Vector{Symbol}()
        interactions = Interactions()
        for (i, bam) in enumerate(bams[r])
            replicate_id = Symbol("$(condition)_$i")
            push!(replicate_ids, replicate_id)
            println("Reading $bam")
            alignments = Alignments(bam; include_secondary_alignments=include_secondary_alignments,
                                    include_alternative_alignments=include_alternative_alignments,
                                    is_reverse_complement=is_reverse_complement)
            println("Annotating alignments...")
            annotate!(alignments, features; prioritize_type=prioritize_type, min_prioritize_overlap=min_prioritize_overlap,
                                            overwrite_type=overwrite_type)
            println("Building graph for replicate $replicate_id...")
            append!(interactions, alignments, replicate_id; min_distance=min_distance, max_ligation_distance=max_ligation_distance,
                filter_types=filter_types, multi_detection_method=multi_detection_method)
            empty!(alignments)
        end
        length(interactions) == 0 && (println("No interactions found!"); continue)

        println("Computing significance levels...")
        addpvalues!(interactions; method=model, include_singles=include_singles, include_read_identity=include_read_identity)
        addpositions!(interactions, features)

        total_reads = sum(interactions.edges[!, :nb_ints])
        above_min_reads = sum(interactions.edges[interactions.edges.nb_ints .>= min_reads, :nb_ints])
        total_ints = nrow(interactions.edges)
        above_min_ints = sum(interactions.edges.nb_ints .>= min_reads)
        total_sig_reads = sum(interactions.edges[interactions.edges.fdr .<= max_fdr, :nb_ints])
        above_min_sig_reads = sum(interactions.edges[(interactions.edges.fdr .<= max_fdr) .& (interactions.edges.nb_ints .>= min_reads), :nb_ints])
        total_sig_ints = sum(interactions.edges.fdr .<= max_fdr)
        above_min_sig_ints = sum((interactions.edges.fdr .<= max_fdr) .& (interactions.edges.nb_ints .>= min_reads))
        println("\n\t\ttotal\treads>=$min_reads\tfdr<=$max_fdr\tboth")
        println("reads:\t\t$total_reads\t$above_min_reads\t\t$total_sig_reads\t\t$above_min_sig_reads")
        println("interactions:\t$total_ints\t$above_min_ints\t\t$total_sig_ints\t\t$above_min_sig_ints\n")
        CSV.write(joinpath(results_path, "interactions", "$(condition).csv"), asdataframe(interactions; output=:edges, min_reads=min_reads, max_fdr=max_fdr))
        CSV.write(joinpath(results_path, "stats", "$(condition).csv"), asdataframe(interactions; output=:stats, min_reads=min_reads, max_fdr=max_fdr))
        CSV.write(joinpath(results_path, "singles", "$(condition).csv"), asdataframe(interactions; output=:nodes, min_reads=min_reads, max_fdr=max_fdr))
    end
    (!overwrite_existing && isfile(joinpath(results_path, "singles.xlsx")) && isfile(joinpath(results_path, "interactions.xlsx"))) && return
    println("Writing tables...")
    singles = CsvFiles(joinpath(results_path, "singles"))
	ints = CsvFiles(joinpath(results_path, "interactions"))
	write(joinpath(results_path, "singles.xlsx"), singles)
	write(joinpath(results_path, "interactions.xlsx"), ints)
end
chimeric_analysis(features::Features, bams::SingleTypeFiles, results_path::String; conditions=conditionsdict(bams),
    filter_types=["rRNA", "tRNA"], min_distance=1000, prioritize_type="sRNA", min_prioritize_overlap=0.8, overwrite_type="IGR", max_ligation_distance=5,
    is_reverse_complement=true, include_secondary_alignments=true, include_alternative_alignments=false, model=:fisher, min_reads=5, max_fdr=0.05,
    overwrite_existing=false, include_read_identity=true, include_singles=true, multi_detection_method=:annotation) =
chimeric_analysis(features, bams, results_path, conditions;
    filter_types=filter_types, min_distance=min_distance, prioritize_type=prioritize_type, min_prioritize_overlap=min_prioritize_overlap,
    overwrite_type=overwrite_type, max_ligation_distance=max_ligation_distance, is_reverse_complement=is_reverse_complement,
    include_secondary_alignments=include_secondary_alignments, include_alternative_alignments=include_alternative_alignments, model=model,
    min_reads=min_reads, max_fdr=max_fdr, overwrite_existing=overwrite_existing, include_read_identity=include_read_identity,
    include_singles=include_singles, multi_detection_method=multi_detection_method)
