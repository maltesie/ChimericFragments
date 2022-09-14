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
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, FASTX.identifier(record)=>FASTA.seqsize(record))
        end
    end
    return interactions_dfs, singles_dfs, stats_dfs, gene_names, genome_info
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

function cytoscape_elements(df::SubDataFrame)
    edges = [Dict("source"=>row.name1, "target"=>row.name2) for row in eachrow(df)]
    nodes = [Dict("id"=>n, "label"=>n) for n in unique(vcat(df.name1, df.name2))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

function circos_data(df::SubDataFrame)
    [Dict("source"=>Dict("id"=>row.ref1, "start"=>row.meanleft1, "end"=>row.meanleft1+1500),
            "target"=>Dict("id"=>row.ref2, "start"=>row.meanleft2, "end"=>row.meanleft2+1500)) for row in eachrow(df)]
end

function table_data(df::SubDataFrame)
    Dict.(pairs.(eachrow(df[!, ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]])))
end

function functional_annotation(gene_names::Vector{String})
end