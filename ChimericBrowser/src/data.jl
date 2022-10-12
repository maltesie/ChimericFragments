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
    end
    gene_names_types = Dict(dname=>merge(Dict(zip(interact.edges.name1, interact.edges.type1)), Dict(zip(interact.edges.name2, interact.edges.type2))) for (dname,interact) in interactions)
    return interactions, gene_names_types, genome_info
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
function cytoscape_elements(df::SubDataFrame, gene_name_type::Dict{String,String})
    total_ints = sum(df.nb_ints)
    edges = [Dict(
        "data"=>Dict(
            "id"=>row.name1*row.name2,
            "source"=>row.name1,
            "target"=>row.name2,
            "current_total"=>total_ints,
            "current_ratio"=>round(row.nb_ints/total_ints; digits=2),
            "interactions"=>row.nb_ints,
            "left1"=>row.left1,
            "right1"=>row.right1,
            "left2"=>row.left2,
            "right2"=>row.right2,
            "rel_int1"=>isnan(row.rel_int1) ? -1.0 : row.rel_int1,
            "rel_lig1"=>isnan(row.rel_lig1) ? -1.0 : row.rel_lig1,
            "rel_int2"=>isnan(row.rel_int2) ? -1.0 : row.rel_int2,
            "rel_lig2"=>isnan(row.rel_lig2) ? -1.0 : row.rel_lig2
        ),
        "classes"=>"srna_edge") for row in eachrow(df)]
    nodes = [Dict(
        "data"=>Dict(
            "id"=>n,
            "label"=>replace(n, ":"=>"\n"),
            "interactions"=>node_sum(df, n),
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2),
            "nb_partners"=>length(union!(Set(df.name1[df.name2 .=== n]), Set(df.name2[df.name1 .=== n]))),
            "current_total"=>total_ints
        ),
        "classes"=>gene_name_type[n]) for n in Set(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame)
    chords_track([Dict{String,Dict{String,Any}}("source"=>Dict("id"=>row.ref1, "start"=>row.left1, "end"=>row.right1),
        "target"=>Dict("id"=>row.ref2, "start"=>row.left2, "end"=>row.right2)) for row in eachrow(df)])
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end

function functional_annotation(gene_names::Vector{String})
end