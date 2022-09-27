function load_data(results_path::String, genome_file::String)
    interactions_path = realpath(joinpath(results_path, "interactions"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".csv")]
    interactions_dfs = Dict(basename(fname)[1:end-4]=>DataFrame(CSV.File(fname; stringtype=String)) for fname in interactions_files)
    gene_names_types = Dict(dname=>merge(Dict(zip(df.name1, df.type1)), Dict(zip(df.name2, df.type2))) for (dname,df) in interactions_dfs)
    singles_path = realpath(joinpath(results_path, "singles"))
    singles_files = [joinpath(singles_path, fname) for fname in readdir(singles_path) if endswith(fname, ".csv")]
    singles_dfs = Dict(basename(fname)[1:end-4]=>DataFrame(CSV.File(fname; stringtype=String)) for fname in singles_files)
    stats_path = realpath(joinpath(results_path,  "stats"))
    stats_files = [joinpath(stats_path, fname) for fname in readdir(stats_path) if endswith(fname, ".csv")]
    stats_matrices = Dict(basename(fname)[1:end-4]=>
        (Dict(row.name1*row.name2=>i for (i,row) in enumerate(eachrow(interactions_dfs[basename(fname)[1:end-4]]))),
        DataFrame(CSV.File(fname; select=collect(1:204)))) for fname in stats_files)
    genome_info = Pair{String,Int}[]
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, FASTX.identifier(record)=>FASTA.seqsize(record))
        end
    end
    return interactions_dfs, singles_dfs, stats_matrices, gene_names_types, genome_info
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
function cytoscape_elements(df::SubDataFrame, srna_type::String, gene_name_type::Dict{String,String})
    total_ints = sum(df.nb_ints)
    srnainteractionindex = (df.type1 .=== srna_type) .!= (df.type2 .=== srna_type)
    edges = [Dict(
        "data"=>Dict(
            "id"=>row.name1*row.name2,
            "source"=>row.name1,
            "target"=>row.name2,
            "current_total"=>total_ints,
            "current_ratio"=>round(row.nb_ints/total_ints; digits=2),
            "interactions"=>row.nb_ints,
            "relpos"=>round(srnainteractionindex[i] && row.type1 === srna_type ? row.relmean2 : row.relmean1; digits=2),
        ),
        "classes"=>(srnainteractionindex[i] ? "srna_edge" : "other_edge")) for (i, row) in enumerate(eachrow(df))]
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
    chords_track([Dict{String,Dict{String,Any}}("source"=>Dict("id"=>row.ref1, "start"=>row.meanleft1, "end"=>row.meanleft1+1500),
        "target"=>Dict("id"=>row.ref2, "start"=>row.meanleft2, "end"=>row.meanleft2+1500)) for row in eachrow(df)])
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end

function functional_annotation(gene_names::Vector{String})
end