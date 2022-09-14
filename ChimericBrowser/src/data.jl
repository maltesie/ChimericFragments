function load_data(results_path::String, genome_file::String)
    interactions_path = realpath(joinpath(results_path, "interactions"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".csv")]
    interactions_dfs = Dict(basename(fname)[1:end-4]=>DataFrame(CSV.File(fname)) for fname in interactions_files)
    gene_names = Dict(dname=>sort(unique(vcat(df.name1, df.name2))) for (dname,df) in interactions_dfs)
    singles_path = realpath(joinpath(results_path, "singles"))
    singles_files = [joinpath(singles_path, fname) for fname in readdir(singles_path) if endswith(fname, ".csv")]
    singles_dfs = Dict(basename(fname)[1:end-4]=>DataFrame(CSV.File(fname)) for fname in singles_files)
    stats_path = realpath(joinpath(results_path,  "stats"))
    stats_files = [joinpath(stats_path, fname) for fname in readdir(stats_path) if endswith(fname, ".csv")]
    stats_dfs = Dict(basename(fname)[1:end-4]=>DataFrame(CSV.File(fname)) for fname in stats_files)
    genome_info = Pair{String,Int}[]
    Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, indentifier(record)=>length(record.sequence))
        end
    end
    return interactions_dfs, singles_dfs, stats_dfs, gene_names, genome_info
end

nthindex(a::Vector{Bool}, n::Int) = n<sum(a) ? sortperm(a; rev=true)[n] : findlast(a) 
function filtered_view(df::DataFrame, search_strings::Vector{String}, min_reads::Int, max_interactions::Int)
    filtered_index = zeros(Bool, nrow(df))
    min_reads_range = 1:findfirst(x->x<min_reads, df.nb_ints)-1
    search_string_index = zeros(Bool, length(min_reads_range))
    for search_string in search_strings
        search_string_index .|= (df.name1[min_reads_range] .=== search_string) .| (df.name2[min_reads_range] .=== search_string)
    end
    n = nthindex(search_string_index, max_interactions)
    filtered_index[1:n] .= search_string_index[1:n]
    return view(df[filtered_index, !])
end

function cytoscape_elements(df::SubDataFrame)
    edges = [Dict("source"=>row.name1, "target"=>row.name2) for row in df]
    nodes = [Dict("id"=>n, "label"=>n) for n in unique(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame)
    [Dict("source"=>Dict("id"=>row.ref1, "start"=>row.mean1, "end"=>row.mean1+1500),
            "target"=>Dict("id"=>row.ref2, "start"=>row.mean2, "end"=>row.mean2+1500)) for row in df]
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end

function functional_annotation(gene_names::Vector{String})
end