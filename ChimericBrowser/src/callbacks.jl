const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks"),
    Output("summary-container", "children"),
    Output("data-tabs", "value"),
    Output("graph", "selectedNodeData"),
    Output("graph", "selectedEdgeData")
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
    State("data-tabs", "value"),
    State("graph", "selectedNodeData"),
    State("graph", "selectedEdgeData")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, bp_len::Int, param_dict::Vector{Pair{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads,
        max_interactions, max_fdr, max_bp_fdr, search_strings, type_strings, layout_value, ligation, exclusive, dataset, tab_value, selected_node, selected_edge
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)
    any(isnothing(v) for v in (min_reads, max_interactions, max_bp_fdr, max_fdr)) && throw(PreventUpdate())
    df = filtered_dfview(interactions[dataset], my_search_strings, my_type_strings, min_reads, max_interactions, Float64(max_fdr),
        Float64(max_bp_fdr), "ligation" in ligation, "exclusive" in exclusive)
    table_output = table_data(df, interactions[dataset])
    cytoscape_output = cytoscape_elements(df, interactions[dataset], layout_value, bp_len, Float64(max_bp_fdr))
    circos_output = circos_data(df, interactions[dataset])
    summary_output = summary_statistics(df, interactions[dataset], param_dict)
    return table_output, cytoscape_output, circos_output, summary_output, tab_value, selected_node, selected_edge
end

const update_dataset_inputs = [
    Input("dropdown-update-dataset", "value")
]
const update_dataset_outputs = [
    Output("min-reads", "value"),
    Output("gene-multi-select", "options"),
    Output("type-multi-select", "options"),
    Output("fdr-value", "value"),
]
const update_dataset_states = [
    State("fdr-value", "value")
]
update_dataset_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, min_reads::Int) =
callback!(app, update_dataset_outputs, update_dataset_inputs, update_dataset_states; prevent_initial_call=false) do dataset, fdr
    return min_reads,
    [Dict("label"=>k, "value"=>k) for k in sort(interactions[dataset].nodes.name)],
    [Dict("label"=>k, "value"=>k) for k in unique(interactions[dataset].nodes.type)],
    fdr
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
const table_column_names = [:name1, :type1, :name2, :type2, :nb_ints, :nb_multi, :in_libs,
                            :fisher_pvalue, :fisher_fdr, :odds_ratio, :bp_pvalue, :bp_fdr, :meanlen1, :nms1, :meanlen2, :nms2]
click_table_button_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}) =
callback!(app, click_table_button_outputs, click_table_button_inputs, click_table_button_states; prevent_initial_call=true) do clicks, dataset,
        min_reads, max_interactions, max_fisher_fdr, max_bp_fdr, search_strings, type_strings, ligation, exclusive
    if clicks>0
        interact = interactions[dataset]
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)
        df = DataFrame(filtered_dfview(interact, my_search_strings, my_type_strings, min_reads, max_interactions,
            Float64(max_fisher_fdr), Float64(max_bp_fdr), "ligation" in ligation, "exclusive" in exclusive))
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

const click_clipboard_inputs = [
    Input("plotly-graph", "clickData")
]
const click_clipboard_outputs = [
    Output("clip-text", "children"),
    Output("clip", "content"),
]
const click_clipboard_states = [
    State("graph", "selectedNodeData"),
    State("graph", "selectedEdgeData"),
]
click_clipboard_callback!(app::Dash.DashApp) =
callback!(app, click_clipboard_outputs, click_clipboard_inputs, click_clipboard_states; prevent_initial_call=true) do click_data, node_data, edge_data
    if !isnothing(click_data) && !isempty(click_data)
        if !isnothing(node_data) && !isempty(node_data)
            p = click_data["points"][1]["x"]
            select_key = click_data["points"][1]["curveNumber"] == 0 ? "lig_as_rna1" : "lig_as_rna2"
            partners = join(["$n: $c" for (n,c) in sort(collect(node_data[1][select_key][p]), by=x->x[2] isa String ? parse(Int, x[2]) : x[2], rev=true)], ", ")
            "<- copy list of partners at position $p", "$partners"
        elseif !isnothing(edge_data) && !isempty(edge_data)
            "<- copy selected basepairing prediction", "$(replace(click_data["points"][1]["text"], "<br>"=>"\n"))"
        end
    else
        throw(PreventUpdate())
    end
end