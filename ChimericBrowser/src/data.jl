struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

function load_data(results_path::String, genome_file::String, min_reads::Int, max_fdr::Float64, max_bp_fdr::Float64)
    interactions_path = realpath(joinpath(results_path, "stats"))
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
        filter!(row -> (row.nb_ints >= min_reads) & (row.fdr <= max_fdr) & (isnan(row.pred_fdr) || (row.pred_fdr .<= max_bp_fdr)), interact.edges)
        interact.edges[!, :meanlen1] = Int.(round.(interact.edges[!, :meanlen1]))
        interact.edges[!, :meanlen2] = Int.(round.(interact.edges[!, :meanlen2]))
        interact.edges[!, :nms1] = round.(interact.edges[!, :nms1], digits=4)
        interact.edges[!, :nms2] = round.(interact.edges[!, :nms2], digits=4)
        interact.edges[:, :name1] = interact.nodes[interact.edges[!,:src], :name]
        interact.edges[:, :name2] = interact.nodes[interact.edges[!,:dst], :name]
        interact.edges[:, :ref1] = interact.nodes[interact.edges[!,:src], :ref]
        interact.edges[:, :ref2] = interact.nodes[interact.edges[!,:dst], :ref]
        interact.edges[:, :type1] = interact.nodes[interact.edges[!,:src], :type]
        interact.edges[:, :type2] = interact.nodes[interact.edges[!,:dst], :type]
        interact.edges[:, :strand1] = interact.nodes[interact.edges[!,:src], :strand]
        interact.edges[:, :strand2] = interact.nodes[interact.edges[!,:dst], :strand]
        interact.edges[:, :in_libs] = sum(eachcol(interact.edges[!, interact.replicate_ids] .!= 0))
        interact.nodes[:, :x] = (rand(rng, nrow(interact.nodes)).+0.5) .* 1200
        interact.nodes[:, :y] = (rand(rng, nrow(interact.nodes)).+0.5) .* 800
        interact.nodes[:, :nb_significant_ints] = zeros(Int, nrow(interact.nodes))
        interact.nodes[:, :nb_significant_partners] = zeros(Int, nrow(interact.nodes))
        for (i,row) in enumerate(eachrow(interact.nodes))
            row[:nb_significant_ints] = sum(interact.edges[(interact.edges.src .== i) .| (interact.edges.dst .== i), :nb_ints])
            names = Set(hcat(interact.nodes[interact.edges[!,:src], :name],interact.nodes[interact.edges[!,:src], :name])[hcat((interact.edges.src .== i),(interact.edges.dst .== i))])
            row[:nb_significant_partners] = length(names)
        end
        sort!(interact.edges, :nb_ints; rev=true)
    end
    gene_name_info = Dict(dname=>Dict(n=>(t,rr,l,r,s,nsingle,nint) for (n,t,l,r,rr,s,nsingle,nint) in
        eachrow(interact.nodes[!, [:name, :type, :left, :right, :ref, :strand, :nb_single, :nb_significant_ints]])) for (dname, interact) in interactions)
    gene_name_position = Dict(dname=>Dict(n=>Dict("x"=>x, "y"=>y) for (n,x,y) in eachrow(interact.nodes[!, [:name, :x, :y]])) for (dname, interact) in interactions)
    return interactions, gene_name_info, gene_name_position, genome_info, genome
end

nthindex(a::BitVector, n::Int) = sum(a)>n ? findall(a)[n] : findlast(a)
function filtered_dfview(df::DataFrame, search_strings::Vector{String}, min_reads::Int, max_interactions::Int,
                            max_fdr::Union{Float64,Int}, max_bp_fdr::Union{Float64,Int}, ligation::Bool, exclusive::Bool)
    filtered_index = falses(nrow(df))
    first_below_min_reads = findfirst(x->x<min_reads, df.nb_ints)
    min_reads_range = 1:(isnothing(first_below_min_reads) ? nrow(df) : first_below_min_reads-1)
    search_string_index = if isempty(search_strings)
        trues(length(min_reads_range))
    else
        s1, s2 = falses(length(min_reads_range)), falses(length(min_reads_range))
        for search_string in search_strings
            s1 .|= (df.name1[min_reads_range] .=== search_string)
            s2 .|= (df.name2[min_reads_range] .=== search_string)
        end
        (exclusive && (length(search_strings) > 1)) ? (s1 .& s2) : (s1 .| s2)
    end
    search_string_index .&= ((df.fdr[min_reads_range] .<= max_fdr) .&
        ((df.pred_fdr[min_reads_range] .<= max_bp_fdr) .| (ligation ? isnan.(df.pred_fdr[min_reads_range]) : falses(length(min_reads_range)))))
    n = nthindex(search_string_index, max_interactions)
    isnothing(n) || (filtered_index[1:n] .= search_string_index[1:n])
    return @view df[filtered_index, :]
end

function get_positions(g::SimpleGraph; xmax=2500, ymax=10000, mean_distance=100, scaling_factor=0.1)
    pos = Vector{Point2}(undef, nv(g))
    components = sort([component for component in connected_components(g)], by=length, rev=true)
    packed_rectangles = GuillotinePacker(xmax, ymax)
    for component in components
        adj = adjacency_matrix(g[component])
        pos[component] .= (stress(adj) .* mean_distance * (1.0 + scaling_factor * sqrt(length(component))))
        minx = minimum(first(p) for p in pos[component])
        maxx = maximum(first(p) for p in pos[component])
        if ((maxx-minx) > (xmax-mean_distance))
            pos[component] .*= (xmax-mean_distance)/(maxx-minx)
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
function clustered_positions(df::SubDataFrame)
    names = collect(Set(vcat(df.name1, df.name2)))
    name_trans = Dict(n=>i for (i,n) in enumerate(names))
    edges = Edge.((name_trans[n1], name_trans[n2]) for (n1, n2) in zip(df.name1, df.name2))
    g = Graph(edges)
    pos = get_positions(g)
    return Dict(names[i]=>Dict("x"=>x, "y"=>y) for (i, (x,y)) in enumerate(pos))
end
node_index(df::SubDataFrame, node_name::String) = (df.name1 .=== node_name) .| (df.name2 .=== node_name)
node_sum(df::SubDataFrame, node_name::String) = sum(df.nb_ints[node_index(df, node_name)])
count_values_rna1(df::SubDataFrame, node_name::String, column_name::Symbol) = collect(counter(df[df.name1 .=== node_name, column_name]))
count_values_rna2(df::SubDataFrame, node_name::String, column_name::Symbol) = collect(counter(df[df.name2 .=== node_name, column_name]))
joint_targets_rna1(df::SubDataFrame, node_name::String, mode_lig::Float64; delim=", ") = join(df[(df.name1 .=== node_name) .& (df.modelig1 .=== mode_lig), :name2], delim)
joint_targets_rna2(df::SubDataFrame, node_name::String, mode_lig::Float64; delim=", ") = join(df[(df.name2 .=== node_name) .& (df.modelig2 .=== mode_lig), :name1], delim)
mode_counts(src::Int, dst::Int, stats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}, mode::Int, dict_index::Int) =
    stats[(src, dst)][dict_index][mode], (mode-1 in keys(stats[(src, dst)][dict_index]) ? stats[(src, dst)][dict_index][mode-1] : 0) + (mode+1 in keys(stats[(src, dst)][dict_index]) ? stats[(src, dst)][dict_index][mode+1] : 0)
function cytoscape_elements(df::SubDataFrame, interact::Interactions, gene_name_info::Dict{String, Tuple{String, String, Int, Int, Char, Int, Int}},
                            gene_name_position::Dict{String, Dict{String, Float64}}, srna_type::String, layout_value::String)
    isempty(df) && return Dict("edges"=>Dict{String,Any}[], "nodes"=>Dict{String,Any}[])
    total_ints = sum(df.nb_ints)
    max_ints = maximum(df.nb_ints)
    srnaindex = hcat(df.type1 .=== srna_type, df.type2 .=== srna_type)
    pos = layout_value == "random" ? gene_name_position : clustered_positions(df)
    edges = [Dict(
        "data"=>Dict(
            "id"=>row.name1*row.name2,
            "source"=>row.name1,
            "target"=>row.name2,
            "current_total"=>total_ints,
            "current_ratio"=>round(row.nb_ints/max_ints; digits=2),
            "interactions"=>row.nb_ints,
            "strand1"=>row.strand1,
            "strand2"=>row.strand2,
            "left1"=>row.left1,
            "right1"=>row.right1,
            "left2"=>row.left2,
            "right2"=>row.right2,
            "ref1"=>row.ref1,
            "ref2"=>row.ref2,
            "len1"=>row.meanlen1,
            "len2"=>row.meanlen2,
            "pred_pvalue"=>(isnan(row.pred_pvalue) ? 2.0 : row.pred_pvalue),
            "modeint1"=>isnan(row.modeint1) ? 0 : row.modeint1,
            "modelig1"=>isnan(row.modelig1) ? 0 : row.modelig1,
            "modeint2"=>isnan(row.modeint2) ? 0 : row.modeint2,
            "modelig2"=>isnan(row.modelig2) ? 0 : row.modelig2,
            "modeintcount1"=>isnan(row.modeint1) ? 0 : interact.edgestats[(row.src, row.dst)][2][row.modeint1],
            "modeligcount1"=>isnan(row.modelig1) ? 0 : interact.edgestats[(row.src, row.dst)][3][row.modelig1],
            "modeintcount2"=>isnan(row.modeint2) ? 0 : interact.edgestats[(row.src, row.dst)][4][row.modeint2],
            "modeligcount2"=>isnan(row.modelig2) ? 0 : interact.edgestats[(row.src, row.dst)][5][row.modelig2],
            "ligcount"=>isnan(row.modelig1) ? 0 : sum(values(interact.edgestats[(row.src, row.dst)][3])),
            "relpos"=> s2 ? (isnan(row.rel_lig1) ? row.rel_int1 : row.rel_lig1) : (isnan(row.rel_lig2) ? row.rel_int2 : row.rel_lig2)
        ),
        "classes"=>s1 != s2 ? "srna_edge" : "other_edge") for (row, (s1, s2)) in zip(eachrow(df), eachrow(srnaindex))]
    nodes = [Dict(
        "data"=>Dict(
            "id"=>n,
            "label"=>replace(n, ":"=>"\n"),
            "interactions"=>node_sum(df, n),
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2),
            "nb_partners"=>length(union!(Set(df.name1[df.name2 .=== n]), Set(df.name2[df.name1 .=== n]))),
            "current_total"=>total_ints,
            "lig_as_rna1"=>Dict("$(Int(gene_name_info[n][5] === '-' ?  gene_name_info[n][4] - k + 1 : k - gene_name_info[n][3] + 1))"=>(
                v, joint_targets_rna1(df, n, k)) for (k,v) in count_values_rna1(df, n, :modelig1) if !isnan(k)),
            "lig_as_rna2"=>Dict("$(Int(gene_name_info[n][5] === '-' ?  gene_name_info[n][4] - k + 1 : k - gene_name_info[n][3] + 1))"=>(
                v, joint_targets_rna2(df, n, k)) for (k,v) in count_values_rna2(df, n, :modelig2) if !isnan(k)),
            "left"=>gene_name_info[n][3],
            "right"=>gene_name_info[n][4],
            "ref"=>gene_name_info[n][2],
            "strand"=>gene_name_info[n][5],
            "nb_single"=>gene_name_info[n][6],
            "nb_ints_total"=>gene_name_info[n][7],
        ),
        "classes"=>gene_name_info[n][1],
        "position"=>pos[n]) for n in Set(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame; min_thickness=1000, max_thickness=3000)
    current_total = sum(df.nb_ints)
    tracks = [chords_track([Dict{String,Any}(
        "source"=>Dict("id"=>row.ref1,
            "start"=>(row.left1 + row.right1)/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(row.left1 + row.right1)/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "target"=>Dict("id"=>row.ref2,
            "start"=>(row.left2 + row.right2)/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(row.left2 + row.right2)/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "data"=>Dict(
                "id"=>row.name1*row.name2,
                "source"=>row.name1,
                "target"=>row.name2,
                "current_total"=>current_total,
                "current_ratio"=>round(row.nb_ints/current_total; digits=2),
                "interactions"=>row.nb_ints,
                "strand1"=>row.strand1,
                "strand2"=>row.strand2,
                "left1"=>row.left1,
                "right1"=>row.right1,
                "left2"=>row.left2,
                "right2"=>row.right2,
                "ref1"=>row.ref1,
                "ref2"=>row.ref2,
                "modeint1"=>isnan(row.modeint1) ? 0 : row.modeint1,
                "modelig1"=>isnan(row.modelig1) ? 0 : row.modelig1,
                "modeint2"=>isnan(row.modeint2) ? 0 : row.modeint2,
                "modelig2"=>isnan(row.modelig2) ? 0 : row.modelig2,
        )
    ) for row in eachrow(df)])]
    return tracks
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end

function interaction_table(df::Union{DataFrame, SubDataFrame}, types::Vector{String})
    types_counter = zeros(Int, (length(types), length(types)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))
    for (t1, t2) in zip(df.type1, df.type2)
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
function summary_statistics(interactions::Interactions, df::SubDataFrame, param_dict::Vector{Pair{String, String}})
    types = sort(unique(interactions.nodes.type))
    dataset_types_table = interaction_table(interactions.edges, types)
    selection_types_table = interaction_table(df, types)
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