struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int,Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64}}
    multichimeras::Dict{Vector{Int}, Int}
    replicate_ids::Vector{Symbol}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

function load_data(results_path::String, genome_file::String, min_reads::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64)
    interactions_path = realpath(joinpath(results_path, "jld"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".jld2")]
    interactions = Dict(basename(fname)[1:end-5]=>Interactions(fname) for fname in interactions_files)
    genome_info = Pair{String,Int}[]
    genome = Dict{String, BioSequences.LongDNA{4}}()
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, identifier(record)=>length(FASTX.sequence(record)))
            push!(genome, identifier(record)=>BioSequences.LongDNA{4}(FASTX.sequence(record)))
        end
    end
    for interact in values(interactions)
        filter!(row -> (row.nb_ints >= min_reads) & (row.fisher_fdr <= max_fisher_fdr) & (isnan(row.bp_fdr) || (row.bp_fdr .<= max_bp_fdr)), interact.edges)
        interact.edges[!, :meanlen1] = Int.(round.(interact.edges[!, :meanlen1]))
        interact.edges[!, :meanlen2] = Int.(round.(interact.edges[!, :meanlen2]))
        interact.edges[!, :nms1] = round.(interact.edges[!, :nms1], digits=4)
        interact.edges[!, :nms2] = round.(interact.edges[!, :nms2], digits=4)
        interact.edges[:, :in_libs] = sum(eachcol(interact.edges[!, interact.replicate_ids] .!= 0))
        interact.nodes[:, :nb_significant_ints] = zeros(Int, nrow(interact.nodes))
        replace!(interact.edges.odds_ratio, Inf=>-1.0)
        for (i,row) in enumerate(eachrow(interact.nodes))
            row[:nb_significant_ints] = sum(interact.edges[(interact.edges.src .== i) .| (interact.edges.dst .== i), :nb_ints])
        end
        sort!(interact.edges, :nb_ints; rev=true)
    end
    return interactions, genome_info, genome
end

nthindex(a::BitVector, n::Int) = sum(a)>n ? findall(a)[n] : findlast(a)
function filtered_dfview(interactions::Interactions, search_strings::Vector{String}, type_strings::Vector{String}, min_reads::Int,
                            max_interactions::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64, ligation::Bool, exclusive::Bool)
    filtered_index = falses(nrow(interactions.edges))
    first_below_min_reads = findfirst(x->x<min_reads, interactions.edges.nb_ints)
    min_reads_range = 1:(isnothing(first_below_min_reads) ? nrow(interactions.edges) : first_below_min_reads-1)
    search_string_index = if isempty(search_strings)
        trues(length(min_reads_range))
    else
        s1, s2 = falses(length(min_reads_range)), falses(length(min_reads_range))
        for search_string in search_strings
            for search_id in findall(interactions.nodes.name .=== search_string)
                s1 .|= (interactions.edges.src[min_reads_range] .=== search_id)
                s2 .|= (interactions.edges.dst[min_reads_range] .=== search_id)
            end
        end
        (exclusive && (length(search_strings) > 1)) ? (s1 .& s2) : (s1 .| s2)
    end
    type_string_index = if isempty(type_strings)
        trues(length(min_reads_range))
    else
        s1, s2 = falses(length(min_reads_range)), falses(length(min_reads_range))
        for type_string in type_strings
            for type_id in findall(interactions.nodes.type .=== type_string)
                s1 .|= (interactions.edges.src[min_reads_range] .=== type_id)
                s2 .|= (interactions.edges.dst[min_reads_range] .=== type_id)
            end
        end
        exclusive ? (s1 .& s2) : (s1 .| s2)
    end
    search_string_index .&= type_string_index
    search_string_index .&= (
        (interactions.edges.fisher_fdr[min_reads_range] .<= max_fisher_fdr) .&
        (
            (interactions.edges.bp_fdr[min_reads_range] .<= max_bp_fdr) .|
            (ligation ? isnan.(interactions.edges.bp_fdr[min_reads_range]) : falses(length(min_reads_range)))
        )
    )
    n = nthindex(search_string_index, max_interactions)
    isnothing(n) || (filtered_index[1:n] .= search_string_index[1:n])
    return @view interactions.edges[filtered_index, :]
end

function make_graph(df::SubDataFrame)
    names = collect(union(Set(df.src), Set(df.dst)))
    name_trans = Dict(n=>i for (i,n) in enumerate(names))
    edges = Edge.((name_trans[n1], name_trans[n2]) for (n1, n2) in zip(df.src, df.dst))
    return Graph(edges), names
end
function get_clustered_positions(g::SimpleGraph; xmax=2500, ymax=10000, mean_distance=100, scaling_factor=0.1)
    pos = Vector{Point2}(undef, nv(g))
    components = sort([component for component in connected_components(g)], by=length, rev=true)
    packed_rectangles = GuillotinePacker(xmax, ymax)
    for component in components
        adj = adjacency_matrix(g[component])
        pos[component] .= (stress(adj) .* mean_distance * (1.0 + scaling_factor * sqrt(length(component))))
        minx = minimum(first(p) for p in pos[component])
        maxx = maximum(first(p) for p in pos[component])
        if ((maxx-minx) > (xmax-mean_distance))
            pos[component] .*= (xmax-mean_distance)/(maxx-minx+1)
            minx = minimum(first(p) for p in pos[component])
            maxx = maximum(first(p) for p in pos[component])
        end
        miny = minimum(last(p) for p in pos[component])
        maxy = maximum(last(p) for p in pos[component])
        pos[component] .-= Point2(minx-(mean_distance/2), miny-(mean_distance/2))
        push!(packed_rectangles, Rect2(0, 0, Int(ceil(maxx-minx+mean_distance)), Int(ceil(maxy-miny+mean_distance))))
    end
    for (component, rect) in zip(components, packed_rectangles.used_rectangles)
        pos[component] .+= Point2(rect.origin[1], rect.origin[2])
    end
    #maxx = maximum(first(p) for p in pos)
    #((xmax-mean_distance) / maxx) > 1.5 && (pos .*= ((xmax-mean_distance) / maxx))
    return pos
end
function clustered_positions(df::SubDataFrame; xmax=2500, ymax=10000, mean_distance=100, scaling_factor=0.1)
    g, names = make_graph(df)
    pos = get_clustered_positions(g; xmax=xmax, ymax=ymax, mean_distance=mean_distance, scaling_factor=scaling_factor)
    return Dict(names[i]=>Dict("x"=>x, "y"=>y) for (i, (x,y)) in enumerate(pos))
end
function grid_positions(df::SubDataFrame; mean_distance=100.0)
    g, names = make_graph(df)
    pos = squaregrid(adjacency_matrix(g))
    pos .-= Point2(minimum(p1 for (p1,_) in pos)-1.0, minimum(p2 for (_,p2) in pos)-1.0)
    pos .*= mean_distance
    return Dict(names[i]=>Dict("x"=>x, "y"=>y) for (i, (x,y)) in enumerate(pos))
end
node_index(df::SubDataFrame, node_id::Int) = (df.src .=== node_id) .| (df.dst .=== node_id)
node_sum(df::SubDataFrame, node_id::Int) = sum(df.nb_ints[node_index(df, node_id)])
function count_ligation_sites_as1(df::SubDataFrame, node_id::Int, interact::Interactions, bp_len::Int, max_fdr::Float64)
    counts = Dict{Int, Float64}()
    partners = Dict{Int, Vector{Tuple{String,Int}}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for partner in df[df.src .== node_id, :dst]
        ligation_points = interact.edgestats[(node_id, partner)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[p][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[p][2]:interact.bpstats[p][3]
                p1 = p[1] + (isnegative ? 1 : -1) * (bp_len - pl)
                if p1 in keys(counts)
                    counts[p1] += c
                    push!(partners[p1], (n, c))
                else
                    counts[p1] = c
                    partners[p1] = [(n, c)]
                end
            end
        end
    end
    return Dict("counts"=>counts, "partners"=>partners)
end
function count_ligation_sites_as2(df::SubDataFrame, node_id::Int, interact::Interactions, bp_len::Int, max_fdr::Float64)
    counts = Dict{Int, Float64}()
    partners = Dict{Int, Vector{Tuple{String,Int}}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for partner in df[df.dst .== node_id, :src]
        ligation_points = interact.edgestats[(partner, node_id)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[p][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[p][4]:interact.bpstats[p][5]
                p2 = p[2] + (isnegative ? -1 : 1) * (bp_len - pl)
                if p2 in keys(counts)
                    counts[p2] += c
                    push!(partners[p2], (n, c))
                else
                    counts[p2] = c
                    partners[p2] = [(n, c)]
                end

            end
        end
    end
    return Dict("counts"=>counts, "partners"=>partners)
end
function cytoscape_elements(df::SubDataFrame, interact::Interactions, layout_value::String, bp_len::Int, max_fdr::Float64)
    isempty(df) && return Dict("edges"=>Dict{String,Any}[], "nodes"=>Dict{String,Any}[])
    total_ints = sum(df.nb_ints)
    max_ints = maximum(df.nb_ints)
    pos = layout_value == "clustered" ? clustered_positions(df) : grid_positions(df)

    edges = [Dict(
        "data"=>Dict(
            "id"=>hash(row.src, hash(row.dst)),
            "source"=>row.src,
            "target"=>row.dst,
            "current_total"=>total_ints,
            "current_ratio"=>round(row.nb_ints/max_ints; digits=2),
            "interactions"=>row.nb_ints,
        ),
        "classes"=>"other_edge"
    ) for row in eachrow(df)]

    nodes = [Dict(
        "data"=>Dict(
            "id"=>n,
            "name"=>interact.nodes.name[n],
            "label"=>replace(interact.nodes.name[n], ":"=>"\n"),
            "interactions"=>node_sum(df, n),
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2),
            "nb_partners"=>length(union!(Set(df.src[df.dst .=== n]), Set(df.dst[df.src .=== n]))),
            "lig_as_rna1"=>count_ligation_sites_as1(df, n, interact, bp_len, max_fdr),
            "lig_as_rna2"=>count_ligation_sites_as2(df, n, interact, bp_len, max_fdr),
        ),
        "classes"=>interact.nodes.type[n],
        "position"=>pos[n]
    ) for n in union(Set(df.src), Set(df.dst))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame, interact::Interactions; min_thickness=1000, max_thickness=3000)
    current_total = sum(df.nb_ints)
    tracks = [chords_track([Dict{String,Any}(
        "source"=>Dict("id"=>interact.nodes.ref[row.src],
            "start"=>(interact.nodes.left[row.src] + interact.nodes.right[row.src])/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(interact.nodes.left[row.src] + interact.nodes.right[row.src])/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "target"=>Dict("id"=>interact.nodes.ref[row.dst],
            "start"=>(interact.nodes.left[row.dst] + interact.nodes.right[row.dst])/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(interact.nodes.left[row.dst] + interact.nodes.right[row.dst])/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "data"=>Dict(
                "id"=>hash(row.src, hash(row.dst)),
                "source"=>row.src,
                "target"=>row.dst,
                "name1"=>interact.nodes.name[row.src],
                "name2"=>interact.nodes.name[row.dst],
                "current_total"=>current_total,
                "current_ratio"=>round(row.nb_ints/current_total; digits=2),
                "interactions"=>row.nb_ints,
        )
    ) for row in eachrow(df)])]
    return tracks
end

function table_data(df::SubDataFrame, interact::Interactions)
    reduced_df = df[:, ["nb_ints", "fisher_fdr", "odds_ratio", "bp_fdr", "in_libs"]]
    reduced_df[:, :name1] = interact.nodes[df[!,:src], :name]
    reduced_df[:, :name2] = interact.nodes[df[!,:dst], :name]
    reduced_df[:, :type1] = interact.nodes[df[!,:src], :type]
    reduced_df[:, :type2] = interact.nodes[df[!,:dst], :type]
    replace!(reduced_df.bp_fdr, NaN => -1.0)
    reduced_df.fisher_fdr = floor.(reduced_df.fisher_fdr; sigdigits=5)
    reduced_df.bp_fdr = floor.(reduced_df.bp_fdr; sigdigits=5)
    reduced_df.odds_ratio = floor.(reduced_df.odds_ratio; sigdigits=5)
    Dict.(pairs.(eachrow(reduced_df)))
end

function interaction_table(df::Union{DataFrame, SubDataFrame}, interactions::Interactions, types::Vector{String})
    types_counter = zeros(Int, (length(types), length(types)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))
    for (t1, t2) in zip(interactions.nodes.type[df.src], interactions.nodes.type[df.dst])
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    table = html_table(children=vcat(
        [html_tr(vcat([html_td("RNA1\\RNA2")], [html_td(h) for h in types]))],
        [html_tr(vcat([html_td(types[j])],[html_td(types_counter[j, i]) for i in eachindex(types)])) for j in eachindex(types)]
    ))
    return table
end

pairs_to_table(dict::Vector{Pair{String, String}}) = html_table([html_tr([html_td(k), html_td(v)]) for (k,v) in dict])
unique_interactions(df::Union{DataFrame, SubDataFrame}) = length(Set(Set((row.src, row.dst)) for row in eachrow(df)))
function summary_statistics(df::SubDataFrame, interactions::Interactions, param_dict::Vector{Pair{String, String}})
    types = sort(unique(interactions.nodes.type))
    dataset_types_table = interaction_table(interactions.edges, interactions, types)
    selection_types_table = interaction_table(df, interactions, types)
    dataset_info = [
        "total interactions:" => "$(nrow(interactions.edges))",
        "unique interactions:" => "$(unique_interactions(interactions.edges))",
        "total annotations" => "$(nrow(interactions.nodes))",
        "interacting annotations" => "$(sum(interactions.nodes.nb_significant_ints .!= 0))",
    ]
    selection_info = [
        "total interactions:" => "$(nrow(df))",
        "unique interactions:" => "$(unique_interactions(df))",
    ]
    html_div(className="horizontal",
        children=[
            html_div([
                html_h3("selection summary:"),
                pairs_to_table(selection_info),
                html_p("annotation types stats:", style=Dict("padding-top"=>"10px")),
                selection_types_table
            ], style=Dict("padding-right"=>"50px")),

            html_div([
                html_h3("dataset summary:"),
                pairs_to_table(dataset_info),
                html_p("annotation types stats:", style=Dict("padding-top"=>"10px")),
                dataset_types_table
            ], style=Dict("border-left"=>"1px solid #000", "padding-left"=>"20px", "padding-right"=>"50px")),

            html_div([
                html_h3("dataset parameter:"),
                pairs_to_table(param_dict)
            ], style=Dict("border-left"=>"1px solid #000", "padding-left"=>"20px")),
        ]
    )
end