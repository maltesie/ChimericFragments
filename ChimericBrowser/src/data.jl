struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int,Dict{Int,Int},Dict{Int,Int},Dict{Int,Int},Dict{Int,Int}}}
    replicate_ids::Vector{Symbol}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

function load_data(results_path::String, genome_file::String)
    interactions_path = realpath(joinpath(results_path, "stats"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".jld2")]
    interactions = Dict(basename(fname)[1:end-5]=>Interactions(fname) for fname in interactions_files)
    genome_info = Pair{String,Int}[]
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, FASTX.identifier(record)=>FASTA.seqsize(record))
        end
    end
    for interact in values(interactions)
        interact.edges[:, :name1] = interact.nodes[interact.edges[!,:src], :name]
        interact.edges[:, :name2] = interact.nodes[interact.edges[!,:dst], :name]
        interact.edges[:, :ref1] = interact.nodes[interact.edges[!,:src], :ref]
        interact.edges[:, :ref2] = interact.nodes[interact.edges[!,:dst], :ref]
        interact.edges[:, :type1] = interact.nodes[interact.edges[!,:src], :type]
        interact.edges[:, :type2] = interact.nodes[interact.edges[!,:dst], :type]
        interact.edges[:, :strand1] = interact.nodes[interact.edges[!,:src], :strand]
        interact.edges[:, :strand2] = interact.nodes[interact.edges[!,:dst], :strand]
        interact.edges[:, :in_libs] = sum(eachcol(interact.edges[!, interact.replicate_ids] .!= 0))
        interact.nodes[:, :x] = rand(nrow(interact.nodes)) .* 100
        interact.nodes[:, :y] = rand(nrow(interact.nodes)) .* 100
        sort!(interact.edges, :nb_ints; rev=true)
    end
    gene_name_type = Dict(dname=>merge(Dict(zip(interact.edges.name1, interact.edges.type1)), Dict(zip(interact.edges.name2, interact.edges.type2))) for (dname,interact) in interactions)
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

node_sum(df::SubDataFrame, node_name::String) = sum(df.nb_ints[(df.name1 .=== node_name) .| (df.name2 .=== node_name)])
nb_partner(df::SubDataFrame, node_name::String) = Set(df.name1[(df.name1 .=== node_name) .| (df.name2 .=== node_name)])
function cytoscape_elements(df::SubDataFrame, gene_name_type::Dict{String,String}, gene_name_position::Dict{String, Dict{String, Float64}}, srna_type::String)
    total_ints = sum(df.nb_ints)
    srnaindex = hcat(df.type1 .=== srna_type, df.type2 .=== srna_type)
    edges = [Dict(
        "data"=>Dict(
            "id"=>row.name1*row.name2,
            "source"=>row.name1,
            "target"=>row.name2,
            "current_total"=>total_ints,
            "current_ratio"=>round(row.nb_ints/total_ints; digits=2),
            "interactions"=>row.nb_ints,
            "strand1"=>row.strand1,
            "strand2"=>row.strand2,
            "left1"=>row.left1,
            "right1"=>row.right1,
            "left2"=>row.left2,
            "right2"=>row.right2,
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
            "interactions"=>node_sum(df, n),
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2),
            "nb_partners"=>length(union!(Set(df.name1[df.name2 .=== n]), Set(df.name2[df.name1 .=== n]))),
            "current_total"=>total_ints
        ),
        "classes"=>gene_name_type[n],
        "position"=>gene_name_position[n]) for n in Set(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame; min_thickness=2000, max_thickness=6000)
    current_total = sum(df.nb_ints)
    chords_track([Dict{String,Any}(
        "nbints"=>row.nb_ints,
        "RNA1"=>row.name1,
        "RNA2"=>row.name2,
        "source"=>Dict("id"=>row.ref1,
            "start"=>(row.left1 + row.right1)/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(row.left1 + row.right1)/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "target"=>Dict("id"=>row.ref2,
            "start"=>(row.left2 + row.right2)/2 - mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(row.left2 + row.right2)/2 + mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2))
    ) for row in eachrow(df)])
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end