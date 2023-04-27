const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks"),
    Output("summary-container", "children"),
    Output("data-tabs", "value")
]
const update_selection_inputs = [
    Input("min-reads", "value"),
    Input("max-interactions", "value"),
    Input("max-fdr", "value"),
    Input("max-bp-fdr", "value"),
    Input("gene-multi-select", "value"),
    Input("type-multi-select", "value"),
    Input("dropdown-update-layout", "value"),
    Input("ligation", "value"),
    Input("exclusive-search", "value"),
]
const update_selection_states = [
    State("dropdown-update-dataset", "value"),
    State("data-tabs", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, sRNA_type::String, param_dict::Vector{Pair{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads,
        max_interactions, max_fdr, max_bp_fdr, search_strings, type_strings, layout_value, ligation, exclusive, dataset, tab_value
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)
    any(isnothing(v) for v in (min_reads, max_interactions, max_bp_fdr, max_fdr)) && throw(PreventUpdate())
    df = filtered_dfview(interactions[dataset], my_search_strings, my_type_strings, min_reads, max_interactions, max_fdr, max_bp_fdr, "ligation" in ligation, "exclusive" in exclusive)
    table_output = table_data(df, interactions[dataset])
    cytoscape_output = cytoscape_elements(df, interactions[dataset], layout_value)
    circos_output = circos_data(df, interactions[dataset])
    summary_output = summary_statistics(df, interactions[dataset], param_dict)
    return table_output, cytoscape_output, circos_output, summary_output, tab_value
end

const update_dataset_inputs = [
    Input("dropdown-update-dataset", "value")
]
const update_dataset_outputs = [
    Output("min-reads", "value"),
    Output("gene-multi-select", "options"),
    Output("type-multi-select", "options"),
]
update_dataset_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, min_reads::Int) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return min_reads,
    [Dict("label"=>k, "value"=>k) for k in sort(interactions[dataset].nodes.name)],
    [Dict("label"=>k, "value"=>k) for k in unique(interactions[dataset].nodes.type)]
end

const update_plots_inputs = [
    Input("fdr-source", "value"),
    Input("fdr-value", "value"),
]
const update_plots_outputs = [
    Output("plot2", "figure"),
    Output("plot1", "figure"),
    Output("plot3", "figure"),
]
const update_plots_states = [
    State("dropdown-update-dataset", "value"),
]
update_plots_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, randseq_model_ecdf::ECDF, genome_model_ecdf::ECDF) =
callback!(app, update_plots_outputs, update_plots_inputs, update_plots_states; prevent_initial_call=false) do plot_type, plot_fdr, dataset
    p1, p2 = plot_pair(interactions[dataset], plot_type, Float64(plot_fdr))
    return bp_score_dist_plot(interactions[dataset], randseq_model_ecdf, genome_model_ecdf, Float64(plot_fdr)), p1, p2
end

normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=0, to_max=100) = Int(floor(to_min + value * (to_max-to_min)))
function cdsframestring(p::Int, idx::Int, interact::Interactions)
    cds, left, right = interact.nodes[idx, [:cds, :left, :right]]
    tp = interact.nodes.strand[idx] == '-' ? (cds > 0 ? cds : right)-p+1 : p-(cds > 0 ? cds : left)+1
    tp <= 0 && (tp -= 1)
    return tp > 0 ? "+$tp" : "$tp"
end

function edge_figure(edge_data::Dash.JSON3.Object, interact::Interactions, genome::Dict{String,BioSequences.LongDNA{4}}, check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel)
    src, dst = parse(Int, edge_data["source"]), parse(Int, edge_data["target"])
    name1, name2 = interact.nodes[src, :name], interact.nodes[dst, :name]
    points1, points2 = [first(p) for p in keys(interact.edgestats[(src,dst)][3])], [last(p) for p in keys(interact.edgestats[(src,dst)][3])]
    tickpos1, tickpos2 = [minimum(points1), maximum(points1)], [minimum(points2), maximum(points2)]
    ticks1, ticks2 = [cdsframestring(p, src, interact) for p in tickpos1], [cdsframestring(p, dst, interact) for p in tickpos2]
    maxints = maximum(values(interact.edgestats[(src, dst)][3]))
    sizes = [ceil(interact.edgestats[(src, dst)][3][p]/maxints*4)*5 for p in zip(points1, points2)]
    bp_plots = [alignment_ascii_plot(src,dst,p1,p2,interact,genome,check_interaction_distances, model) for (p1,p2) in zip(points1, points2)]
    colors = adjust(PValues([interact.bpstats[p][1] for p in zip(points1, points2)]), BenjaminiHochberg())
    return plot(scatter(y=points1, x=points2, mode="markers", marker=attr(size=sizes, color=colors, colorbar=attr(title="FDR", orientation="h"),
            colorscale="Reds", reversescale=true, cmin=0.0, cmax=1.0), name = "ligation points",
            text=bp_plots, hoverinfo="text", hovertemplate="<span style=\"font-family:'Lucida Console', monospace\">%{text}</span><extra></extra>"),
        Layout(title = "$(edge_data["interactions"]) interactions",
            yaxis=attr(title="RNA1: $name1", tickmode="array", tickvals=tickpos1, ticktext=ticks1),
            xaxis=attr(title="RNA2: $name2", tickmode="array", tickvals=tickpos2, ticktext=ticks2)))
end

function node_figure(node_data::Dash.JSON3.Object, interact::Interactions)
    idx = parse(Int, node_data["id"])
    name = interact.nodes.name[idx]
    minpos = minimum(minimum(parse(Int, String(k)) for (k,_) in node_data[select_key]) for select_key in ("lig_as_rna1", "lig_as_rna2") if length(node_data[select_key])>0)
    maxpos = maximum(maximum(parse(Int, String(k)) for (k,_) in node_data[select_key]) for select_key in ("lig_as_rna1", "lig_as_rna2") if length(node_data[select_key])>0)
    tickpos = [minpos, maxpos]
    ticktext = [cdsframestring(p, idx, interact) for p in tickpos]
    return plot([
            begin
                ligationpoints = [parse(Int, String(k))=>v for (k,v) in node_data[select_key]]
                kv = sort(ligationpoints, by=x->x[1])
                positions, counts = Int[t[1] for t in kv], Int[t[2] for t in kv]

                ticks = [cdsframestring(p, idx, interact) for p in tickpos]
                scatter(x = positions, y = counts, fill="tozeroy", name = legend)
            end
            for (select_key, legend) in zip(("lig_as_rna1", "lig_as_rna2"), ("as RNA1", "as RNA2"))
        ],
        Layout(title = "$name ($(node_data["nb_partners"]) partners)", xaxis_title = "position", yaxis_title = "count",
            xaxis=attr(tickmode="array", tickvals=tickpos, ticktext=ticktext), showlegend=false)
    )
end

const empty_figure = (
    data = [
        (x = [], y = [], type = "scatter", name = ""),
    ],
    layout = (title = "Please select an edge<br>or a node in the graph.",)
)

const update_selected_element_inputs = [
    Input("graph", "selectedNodeData"),
    Input("graph", "selectedEdgeData")
]
const update_selected_element_outputs = [
    Output("plotly-graph", "figure"),
]
const update_selected_element_states = [
    State("dropdown-update-dataset", "value")
]
update_selected_element_callback!(app::Dash.DashApp, genome::Dict{String,BioSequences.LongDNA{4}}, interactions::Dict{String, Interactions},
        check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs, update_selected_element_states; prevent_initial_call=true) do node_data, edge_data, dataset

    no_node_data = isnothing(node_data) || isempty(node_data)
    no_edge_data = isnothing(edge_data) || isempty(edge_data)
    no_node_data && no_edge_data && return [empty_figure]
    no_edge_data && return [node_figure(node_data[1], interactions[dataset])]
    return [edge_figure(edge_data[1], interactions[dataset], genome, check_interaction_distances, model)]
end

click_cyto_button_callback!(app::Dash.DashApp) =
callback!(app, Output("graph", "generateImage"), Input("save-svg", "n_clicks"), State("dropdown-update-dataset", "value"); prevent_initial_call=true) do clicks, dataset
    clicks>0 && return Dict("type"=>"svg", "action"=>"download", "filename"=>"$(dataset)_graph")
end

const click_table_button_inputs = [
    Input("btn-csv", "n_clicks")
]
const click_table_button_outputs = [
    Output("download-dataframe-csv", "data")
]
const click_table_button_states = [
    State("dropdown-update-dataset", "value"),
    State("min-reads", "value"),
    State("max-interactions", "value"),
    State("max-fdr", "value"),
    State("max-bp-fdr", "value"),
    State("gene-multi-select", "value"),
    State("type-multi-select", "value"),
    State("ligation", "value"),
    State("exclusive-search", "value"),
]
const table_column_names = [:name1, :type1, :name2, :type2, :nb_ints, :nb_multi, :in_libs, :pvalue, :fdr, :pred_pvalue, :pred_fdr,
                            :modeint1, :modelig1, :meanlen1, :nms1, :modeint2, :modelig2, :meanlen2, :nms2]
click_table_button_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}) =
callback!(app, click_table_button_outputs, click_table_button_inputs, click_table_button_states; prevent_initial_call=true) do clicks, dataset,
        min_reads, max_interactions, max_fdr, max_bp_fdr, search_strings, type_strings, ligation, exclusive
    if clicks>0
        interact = interactions[dataset]
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)
        df = DataFrame(filtered_dfview(interact, my_search_strings, my_type_strings, min_reads, max_interactions, max_fdr, max_bp_fdr,"ligation" in ligation, "exclusive" in exclusive))
        df.name1 = interact.nodes.name[df.src]
        df.name2 = interact.nodes.name[df.dst]
        df.type1 = interact.nodes.type[df.src]
        df.type2 = interact.nodes.type[df.dst]
        df = df[!, table_column_names]
        csvrowwriteriterator = CSV.RowWriter(df)
        dfstring = join(collect(csvrowwriteriterator))
        return [Dict("filename"=>"$(dataset)_table.csv", "content"=>dfstring ,"base64"=>false)]
    end
end

const click_add_node_inputs = [
    Input("add-selected-btn", "n_clicks")
]
const click_add_node_outputs = Output("gene-multi-select", "value")

const click_add_node_states = [
    State("graph", "selectedNodeData"),
    State("gene-multi-select", "value"),
]
click_add_node_callback!(app::Dash.DashApp) =
callback!(app, click_add_node_outputs, click_add_node_inputs, click_add_node_states; prevent_initial_call=true) do clicks, node_data, search_strings
    if clicks>0 && !isnothing(node_data) && !isempty(node_data)
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        node_data[1]["name"] in my_search_strings && throw(PreventUpdate())
        return push!(my_search_strings, node_data[1]["name"])
    else
        throw(PreventUpdate())
    end
end