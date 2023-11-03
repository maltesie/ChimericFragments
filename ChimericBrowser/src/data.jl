# data.jl defines functions to load data from JLD2 files, filter the data and generate objects compatible with the different
# visualization modes

# Interactions struct is the data container created and populated by analyze.jl and ChimericAnalysis. It holds all relevant information
# on interactions, genes and ligation points
struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int, Int, Int, Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int64}}
    multichimeras::Dict{Vector{Int}, Int}
    replicate_ids::Vector{Symbol}
    counts::Dict{Symbol,Vector{Int}}
end

# Helper method to create an Interactions struct from a .jld2 file
Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

# Initial function to prepare the data from .jld2 files for use in the browser.
function load_data(results_path::String, genome_file::String, min_reads::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64)
    # look for all .jld2 files in the results folder
    interactions_path = realpath(joinpath(results_path, "jld"))
    interactions_files = [joinpath(interactions_path, fname) for fname in readdir(interactions_path) if endswith(fname, ".jld2")]
    # load files into Interactions struct and use filenames as labels for the contained datasets
    interactions = Dict(basename(fname)[1:end-5]=>Interactions(fname) for fname in interactions_files)
    # read genome file and collect info for circos plots
    genome_info = Pair{String,Int}[]
    genome = Dict{String, BioSequences.LongDNA{4}}()
    FASTA.Reader(open(genome_file)) do reader
        for record in reader
            push!(genome_info, identifier(record)=>length(FASTX.sequence(record)))
            push!(genome, identifier(record)=>BioSequences.LongDNA{4}(FASTX.sequence(record)))
        end
    end
    # adjust Float precision for use in the browser
    for interact in values(interactions)
        filter!(row -> (row.nb_ints >= min_reads) & (row.fisher_fdr <= max_fisher_fdr) & (isnan(row.bp_fdr) || (row.bp_fdr .<= max_bp_fdr)), interact.edges)
        interact.edges[!, :meanlen1] = Int.(round.(interact.edges[!, :meanlen1]))
        interact.edges[!, :meanlen2] = Int.(round.(interact.edges[!, :meanlen2]))
        interact.edges[!, :nms1] = round.(interact.edges[!, :nms1], digits=4)
        interact.edges[!, :nms2] = round.(interact.edges[!, :nms2], digits=4)
        interact.edges[:, :in_libs] = sum(eachcol(interact.edges[!, interact.replicate_ids] .!= 0))
        interact.nodes[:, :nb_significant_ints] = zeros(Int, nrow(interact.nodes))
        # handle Inf values, Dash uses JSON3 to pass data around and it does not like Infs
        replace!(interact.edges.odds_ratio, Inf=>-1.0)
        # compute number of interactions for a gene depending from the total dataset
        for (i,row) in enumerate(eachrow(interact.nodes))
            row[:nb_significant_ints] = sum(interact.edges[(interact.edges.src .== i) .| (interact.edges.dst .== i), :nb_ints])
        end
        # sort interactions table to work with filtering methods
        sort!(interact.edges, :nb_ints; rev=true)
    end
    return interactions, genome_info, genome
end

# auxiliary function, detemines the index of the first interaction after a certain amount of interactions according to selected criteria is reached
nthindex(a::BitVector, n::Int) = sum(a)>n ? findall(a)[n] : findlast(a)
# main search and filter function, returns a DataFrame view with entries according to selected criteria
function filtered_dfview(interactions::Interactions, search_strings::Vector{String}, type_strings::Vector{String}, min_reads::Int,
                            max_interactions::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64, ligation::Bool, exclusive::Bool)
    #initialize filter index
    filtered_index = falses(nrow(interactions.edges))
    # set cutoff for minimum number of reads per interaction
    first_below_min_reads = findfirst(x->x<min_reads, interactions.edges.nb_ints)
    min_reads_range = 1:(isnothing(first_below_min_reads) ? nrow(interactions.edges) : first_below_min_reads-1)

    # handle searching for specific genes
    search_string_index = if isempty(search_strings)
        trues(length(min_reads_range))
    else
        # search for specified genes in the interactions
        s1, s2 = falses(length(min_reads_range)), falses(length(min_reads_range))
        for search_string in search_strings
            # get internal id of the current gene
            for search_id in findall(interactions.nodes.name .=== search_string)
                # index of all interactions where the gene is on position 1
                s1 .|= (interactions.edges.src[min_reads_range] .=== search_id)
                # index of all interactions where the gene is on position 2
                s2 .|= (interactions.edges.dst[min_reads_range] .=== search_id)
            end
        end
        # search behavior can be adjusted to AND or OR the selected criteria in the control called "exclusive search"
        (exclusive && (length(search_strings) > 1)) ? (s1 .& s2) : (s1 .| s2)
    end

    # do the same as above for annotation types instead of gene names
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
        different_type = (length(type_strings) > 1) ?
            interactions.nodes.type[interactions.edges.src[min_reads_range]] .!== interactions.nodes.type[interactions.edges.dst[min_reads_range]] :
            trues(length(min_reads_range))
        exclusive ? (s1 .& s2 .& different_type) : (s1 .| s2)
    end
    # combine search criteria
    search_string_index .&= type_string_index
    search_string_index .&= (
        (interactions.edges.fisher_fdr[min_reads_range] .<= max_fisher_fdr) .&
        (
            (interactions.edges.bp_fdr[min_reads_range] .<= max_bp_fdr) .|
            (ligation ? isnan.(interactions.edges.bp_fdr[min_reads_range]) : falses(length(min_reads_range)))
        )
    )
    # cut after max number of interactions is reached
    n = max_interactions > 0 ? nthindex(search_string_index, max_interactions) : nothing
    isnothing(n) || (filtered_index[1:n] .= search_string_index[1:n])
    return @view interactions.edges[filtered_index, :]
end

# generate a Graph defined in LightGraphs dependency
function make_graph(df::SubDataFrame)
    names = collect(union(Set(df.src), Set(df.dst)))
    name_trans = Dict(n=>i for (i,n) in enumerate(names))
    edges = Edge.((name_trans[n1], name_trans[n2]) for (n1, n2) in zip(df.src, df.dst))
    return Graph(edges), names
end
# compute positions for placing nodes for a given graph
function get_clustered_positions(g::SimpleGraph; xmax=2500, ymax=10000, mean_distance=100, scaling_factor=0.1)
    pos = Vector{Point2}(undef, nv(g))
    # collect and sort connected components of the graph of interactions
    components = sort([component for component in connected_components(g)], by=length, rev=true)
    # GuillotinePacker finds good placements for rectangles in a larger box with minimal free space between them
    packed_rectangles = GuillotinePacker(xmax, ymax)
    for component in components
        # compute the adjacency matrix of a connected component
        adj = adjacency_matrix(g[component])
        # run stress majorization on the adjacency matrix to generate relative positions within the component
        pos[component] .= (stress(adj) .* mean_distance * (1.0 + scaling_factor * sqrt(length(component))))
        # compute boundary box of the positions of the component
        minx = minimum(first(p) for p in pos[component])
        maxx = maximum(first(p) for p in pos[component])
        # scale the positions of a component down if it is to large
        if ((maxx-minx) > (xmax-mean_distance))
            pos[component] .*= (xmax-mean_distance)/(maxx-minx+1)
            minx = minimum(first(p) for p in pos[component])
            maxx = maximum(first(p) for p in pos[component])
        end
        # recompute boundary box of the positions of the component
        miny = minimum(last(p) for p in pos[component])
        maxy = maximum(last(p) for p in pos[component])
        # add margin around the component
        pos[component] .-= Point2(minx-(mean_distance/2), miny-(mean_distance/2))
        # send to GuillotinePacker to find coordinate on canvas
        push!(packed_rectangles, Rect2(0, 0, Int(ceil(maxx-minx+mean_distance)), Int(ceil(maxy-miny+mean_distance))))
    end
    # collect points from GuillotinePacker
    for (component, rect) in zip(components, packed_rectangles.used_rectangles)
        pos[component] .+= Point2(rect.origin[1], rect.origin[2])
    end
    return pos
end
# handle graph generation and positioning with stress majorization for a filtered set of interactions
function clustered_positions(df::SubDataFrame; xmax=2500, ymax=10000, mean_distance=100, scaling_factor=0.1)
    g, names = make_graph(df)
    pos = get_clustered_positions(g; xmax=xmax, ymax=ymax, mean_distance=mean_distance, scaling_factor=scaling_factor)
    return Dict(names[i]=>Dict("x"=>x, "y"=>y) for (i, (x,y)) in enumerate(pos))
end
# handle graph generation and positioning on a grid for a filtered set of interactions
function grid_positions(df::SubDataFrame; mean_distance=100.0)
    g, names = make_graph(df)
    pos = squaregrid(adjacency_matrix(g))
    pos .-= Point2(minimum(p1 for (p1,_) in pos)-1.0, minimum(p2 for (_,p2) in pos)-1.0)
    pos .*= mean_distance
    return Dict(names[i]=>Dict("x"=>x, "y"=>y) for (i, (x,y)) in enumerate(pos))
end

# auxiliary functions for generating the data compatible with the cytoscape graph drawing component
node_index(df::SubDataFrame, node_id::Int) = (df.src .=== node_id) .| (df.dst .=== node_id)
node_sum(df::SubDataFrame, node_id::Int) = sum(df.nb_ints[node_index(df, node_id)])

# compute all partners with a ligation point with the selected gene if the gene is on position 1 in the pair
lig_as_rna1(df::SubDataFrame, node_id::Int, interact::Interactions) =
    [partner for partner in df[df.src .== node_id, :dst] if length(interact.edgestats[(node_id, partner)][3]) > 0]

# compute all partners with a ligation point with the selected gene if the gene is on position 2 in the pair
lig_as_rna2(df::SubDataFrame, node_id::Int, interact::Interactions) =
    [partner for partner in df[df.dst .== node_id, :src] if length(interact.edgestats[(partner, node_id)][3]) > 0]

#main function for generating the data compatible with the cytoscape graph drawing component
function cytoscape_elements(df::SubDataFrame, interact::Interactions, layout_value::String)

    # handle empty results of a search
    isempty(df) && return Dict("edges"=>Dict{String,Any}[], "nodes"=>Dict{String,Any}[])

    total_ints = sum(df.nb_ints)
    max_ints = maximum(df.nb_ints)

    # compute positions for the graph of currently selected interactions
    pos = layout_value == "clustered" ? clustered_positions(df) : grid_positions(df)

    # edges object with information for directed edges (source->target)
    edges = [Dict(
        "data"=>Dict(
            "id"=>hash(row.src, hash(row.dst)), # internal id of the edge
            "source"=>row.src,  # draw from node with source id
            "target"=>row.dst,  # to node with target id
            "current_total"=>total_ints,    # sum of reads in currently selected interactions
            "current_ratio"=>round(row.nb_ints/max_ints; digits=2), # ratio is used to set thickness of the edge
            "bp_fdr"=>isnan(row.bp_fdr) ? 2.0 : row.bp_fdr, # reformat fdr, Dash cant handle NaN
            "interactions"=>row.nb_ints,
        ),
        "classes"=>isnan(row.bp_fdr) ? "no_ligation_edge" : "ligation_edge"
    ) for row in eachrow(df)]

    # nodes object with local network information for for genes in the current selection of interactions
    nodes = [Dict(
        "data"=>Dict(
            "id"=>n, # internal id of the node, edges use them as source and target
            "name"=>interact.nodes.name[n], # internal name of the node
            "label"=>replace(interact.nodes.name[n], ":"=>"\n"), # reformat labels of intergenic regions to two lines
            "interactions"=>node_sum(df, n), # sum of reads per interaction in the current network
            "current_ratio"=>round(node_sum(df, n)/total_ints; digits=2), # ratio is used to set the size of the node
            "nb_partners"=>length(union!(Set(df.src[df.dst .=== n]), Set(df.dst[df.src .=== n]))), # Number of partners in the current selection
            "lig_as_rna1"=>lig_as_rna1(df, n, interact), # info for plotting aggregated baspairing for this node
            "lig_as_rna2"=>lig_as_rna2(df, n, interact), # info for plotting aggregated baspairing for this node
        ),
        "classes"=>interact.nodes.type[n], # define how the node is drawn
        "position"=>pos[n] # positions of the nodes
    ) for n in union(Set(df.src), Set(df.dst))]
    return Dict("edges"=>edges, "nodes"=>nodes)
end

# main function to generate data compatible with the circos plot component
function circos_data(df::SubDataFrame, interact::Interactions; min_thickness=1000, max_thickness=3000)
    current_total = sum(df.nb_ints)
    tracks = [chords_track([Dict{String,Any}(
        # the circos plot draws undirected lines from source to target
        # start and end define the thickness of the line which is relative to the number of reads the interaction was found in
        "source"=>Dict("id"=>interact.nodes.ref[row.src],
            "start"=>(interact.nodes.left[row.src] + interact.nodes.right[row.src])/2 -
                mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(interact.nodes.left[row.src] + interact.nodes.right[row.src])/2 +
                mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        "target"=>Dict("id"=>interact.nodes.ref[row.dst],
            "start"=>(interact.nodes.left[row.dst] + interact.nodes.right[row.dst])/2 -
                mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2),
            "end"=>(interact.nodes.left[row.dst] + interact.nodes.right[row.dst])/2 +
                mapvalue(row.nb_ints/current_total; to_min=min_thickness/2, to_max=max_thickness/2)),
        # additional data can be used to layout the circos plot. currently not used
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

# main function to generate data compatible with the data table component of Dash
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

# html table for stats of interacting annotation types
function interaction_table(df::Union{DataFrame, SubDataFrame}, interactions::Interactions, types::Vector{String})
    # matrix of size n x n, with n the number of annotation types present in the dataset
    types_counter = zeros(Int, (length(types), length(types)))
    type_trans = Dict{String, Int}(t=>i for (i,t) in enumerate(types))
    # count pairs of annotation types found in interactions
    for (t1, t2) in zip(interactions.nodes.type[df.src], interactions.nodes.type[df.dst])
        types_counter[type_trans[t1], type_trans[t2]] += 1
    end
    # create html table using Dash html elements
    table = html_table(children=vcat(
        [html_tr(vcat([html_td("RNA1\\RNA2")], [html_td(h) for h in types]))],
        [html_tr(vcat([html_td(types[j])],[html_td(types_counter[j, i]) for i in eachindex(types)])) for j in eachindex(types)]
    ))
    return table
end

# html table for summary of annotation types
singles_table(types::Vector{String}, interactions::Interactions) = html_table(
    vcat(
        [html_tr([html_td("type"), html_td("annotations"), html_td("reads")])],
        [html_tr([html_td(t), html_td(sum(interactions.nodes.type .== t)), html_td(sum(interactions.nodes.nb_single[interactions.nodes.type .== t]))]) for t in types]
    )
)

# table of replecate correlation values, currently not in use
const row_names = ["single:", "unclassified:", "selfchimera:", "excluded:", "chimeria:", "multichimera:"]
replicate_table(interactions::Interactions) = html_table(
    vcat(
        [html_tr(vcat([html_td("")], [html_td(String(s)) for s in interactions.replicate_ids]))],
        [
            html_tr(vcat([row_names[i]], [html_td(v) for v in [interactions.counts[replicate_id][i] for replicate_id in interactions.replicate_ids]])) for i in 1:6
        ]
    )
)

# auxiliary functions to generate tables
pairs_to_table(dict::Vector{Pair{String, String}}) = html_table([html_tr([html_td(k), html_td(v)]) for (k,v) in dict])
unique_interactions(df::Union{DataFrame, SubDataFrame}) = length(Set(Set((row.src, row.dst)) for row in eachrow(df)))
function summary_statistics(df::SubDataFrame, interactions::Interactions, param_dict::Vector{Pair{String, String}})
    types = sort(unique(interactions.nodes.type))
    # generate table for total dataset
    dataset_types_table = interaction_table(interactions.edges, interactions, types)
    # generate table for current selection of interactions
    selection_types_table = interaction_table(df, interactions, types)
    # info for total dataset
    dataset_info = [
        "total interactions:" => "$(nrow(interactions.edges))",
        "unique interactions:" => "$(unique_interactions(interactions.edges))",
        "total annotations" => "$(nrow(interactions.nodes))",
        "interacting annotations" => "$(sum(interactions.nodes.nb_significant_ints .!= 0))",
    ]
    # info for currently selected interactions
    selection_info = [
        "total interactions:" => "$(nrow(df))",
        "unique interactions:" => "$(unique_interactions(df))",
        "total annotations" => "$(nrow(interactions.nodes))",
        "interacting annotations" => "$(length(unique(vcat(df.src, df.dst))))",
    ]
    # html object layouting the summary tab containing a tabular view for the total dataset, the current selection of interactions
    # and a summary of important parameters set in the config
    html_div(className="horizontal",
        children=[
            html_div([
                html_h3("selection summary:"),
                html_p("interaction stats:", style=Dict("padding-top"=>"10px", "font-weight"=>"bold")),
                pairs_to_table(selection_info),
                selection_types_table
            ], style=Dict("padding-right"=>"50px")),

            html_div([
                html_h3("dataset summary:"),
                html_p("interaction stats:", style=Dict("padding-top"=>"10px", "font-weight"=>"bold")),
                pairs_to_table(dataset_info),
                dataset_types_table,
                html_p("single stats:", style=Dict("padding-top"=>"10px", "font-weight"=>"bold")),
                singles_table(types, interactions),
            ], style=Dict("border-left"=>"1px solid #000", "padding-left"=>"20px", "padding-right"=>"50px")),

            html_div([
                html_h3("dataset parameter:"),
                pairs_to_table(param_dict)
            ], style=Dict("border-left"=>"1px solid #000", "padding-left"=>"20px")),
        ]
    )
end