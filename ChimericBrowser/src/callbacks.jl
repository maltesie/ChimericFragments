const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks")
]
const update_selection_inputs = [
    Input("reads-slider", "value"),
    Input("max-interactions", "value"),
    Input("gene-multi-select", "value"),
]
const update_selection_states = [
    State("dropdown-update-dataset", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions_dfs::Dict{String,DataFrame}, srna_type::String, gene_name_types::Dict{String, Dict{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, search_strings, dataset
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    df = filtered_dfview(interactions_dfs[dataset], my_search_strings, min_reads, max_interactions)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df, srna_type, gene_name_types[dataset])
    circos_output = circos_data(df)
    return table_output, cytoscape_output, circos_output
end

const update_dataset_inputs = [
    Input("dropdown-update-dataset", "value")
]
const update_dataset_outputs = [
    Output("reads-slider", "value"),
    Output("gene-multi-select", "options")
]
update_dataset_callback!(app::Dash.DashApp, gene_names_types::Dict{String,Dict{String,String}}) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return 0, [Dict("label"=>k, "value"=>k) for k in sort(collect(keys(gene_names_types[dataset])))]
end

function edge_info(stats_df::Tuple{Dict{String,Int}, DataFrame}, source_name::String, target_name::String, interactions::Int)
    stats_row = stats_df[2][stats_df[1][source_name*target_name], :]
    return [html_div(id="edge-info", children=[
        "$source_name and $target_name were found in $interactions chimeric read" * (interactions>1 ? "s" : "") * ".",
        "$source_name fragment distribution:", string(stats_row[1:10]),
        "$target_name fragment distribution:", string(stats_row[end-9:end])
    ])]
end

const update_selected_element_inputs = [
    Input("graph", "selectedNodeData"),
    Input("graph", "selectedEdgeData")
]
const update_selected_element_outputs = [
    Output("info-output", "children")
]
const update_selected_element_states = [
    State("dropdown-update-dataset", "value")
]
update_selected_element_callback!(app::Dash.DashApp, stats_dfs::Dict{String, Tuple{Dict{String,Int}, DataFrame}}) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs, update_selected_element_states; prevent_initial_call=true) do node_data, edge_data, dataset
    no_node_data = isnothing(node_data) || isempty(node_data)
    no_edge_data = isnothing(edge_data) || isempty(edge_data)
    no_node_data && no_edge_data && return ["Select an edge or node in the graph to display additional information."]
    no_edge_data && return ["$(node_data[1]["id"]) has $(node_data[1]["nb_partners"]) partner" * (node_data[1]["nb_partners"]>1 ? "s" : "") * " in the current selection."]
    return edge_info(stats_dfs[dataset], edge_data[1]["source"], edge_data[1]["target"], edge_data[1]["interactions"])
end