struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

function load_data(results_path::String, genome_file::String, min_reads::Int, max_fdr::Float64)
    interactions_path = realpath(joinpath(results_path, "stats"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".jld2")]
    interactions = Dict(basename(fname)[1:end-5]=>Interactions(fname) for fname in interactions_files)
    genome_info = Pair{String,Int}[]
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, identifier(record)=>length(sequence(record)))
        end
    end
    for interact in values(interactions)
        filter!(row -> (row.nb_ints >= min_reads) & (row.fdr <= max_fdr), interact.edges)
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
        sort!(interact.edges, :nb_ints; rev=true)
    end
    gene_name_info = Dict(dname=>Dict(n=>(t,rr,l,r) for (n,t,l,r,rr) in interact.nodes[!, [:name, :type, :left, :right, :ref]]) for (dname,interact) in interactions)
    gene_name_position = Dict(dname=>Dict(n=>Dict("x"=>x, "y"=>y) for (n,x,y) in zip(interact.nodes.name, interact.nodes.x, interact.nodes.y)) for (dname,interact) in interactions)
    return interactions, gene_name_type, gene_name_position, genome_info
end

nthindex(a::Vector{Bool}, n::Int) = sum(a)>n ? findall(a)[n] : findlast(a)
function filtered_dfview(df::DataFrame, search_strings::Vector{String}, min_reads::Int, max_interactions::Int)
    filtered_index = zeros(Bool, nrow(df))
    first_below_min_reads = findfirst(x->x<min_reads, df.nb_ints)
    min_reads_range = 1:(isnothing(first_below_min_reads) ? nrow(df) : first_below_min_reads-1)
    search_string_index = isempty(search_strings) ? ones(Bool, length(min_reads_range)) : zeros(Bool, length(min_reads_range))
    for search_string in search_strings
        search_string_index .|= ((df.name1[min_reads_range] .=== search_string) .| (df.name2[min_reads_range] .=== search_string))
    end
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
function cytoscape_elements(df::SubDataFrame, gene_name_info::Dict{String,String}, gene_name_position::Dict{String, Dict{String, Float64}}, srna_type::String, layout_value::String)
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
            "modeint1"=>isnan(row.modeint1) ? 0 : row.modeint1,
            "modelig1"=>isnan(row.modelig1) ? 0 : row.modelig1,
            "modeint2"=>isnan(row.modeint2) ? 0 : row.modeint2,
            "modelig2"=>isnan(row.modelig2) ? 0 : row.modelig2,
            "relpos"=> s2 ? (isnan(row.rel_lig1) ? row.rel_int1 : row.rel_lig1) : (isnan(row.rel_lig2) ? row.rel_int2 : row.rel_lig2)
        ),
        "classes"=>s1 != s2 ? "srna_edge" : "other_edge") for (row, (s1, s2)) in zip(eachrow(df), eachrow(srnaindex))]
    nodes = [Dict(
        "data"=>Dict(
            "id"=>n,
            "label"=>replace(n, ":"=>"\n"),
            "left"=>
            "interactions"=>node_sum(df, n),
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2),
            "nb_partners"=>length(union!(Set(df.name1[df.name2 .=== n]), Set(df.name2[df.name1 .=== n]))),
            "current_total"=>total_ints,
            "lig_as_rna1"=>Dict("$(Int(k))"=>v for (k,v) in sort(count_values_rna1(df, n, :modelig1), by=x->x[2], rev=true) if !isnan(k)),
            "lig_as_rna2"=>Dict("$(Int(k))"=>v for (k,v) in sort(count_values_rna2(df, n, :modelig2), by=x->x[2], rev=true) if !isnan(k)),
        ),
        "classes"=>gene_name_info[n][1],
        "position"=>pos[n]) for n in Set(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame; min_thickness=2000, max_thickness=5000)
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
