test_callback!(app::Dash.DashApp) =
callback!(app, Output("info-output", "children"), Input("reads-slider", "value")) do slider_value
    "$slider_value"
end

const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks"),
    Output("info-output", "children")
]

const update_selection_inputs = [
    Input("dropdown-update-dataset", "value"),
    Input("reads-slider", "value"),
    Input("max-interactions", "value"),
    Input("gene-multi-select", "value"),
]

update_selection_callback!(app::Dash.DashApp, interactions_dfs::Dict{String,DataFrame}, srna_type::String, gene_name_types::Dict{String, Dict{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs) do dataset, min_reads, max_interactions, search_strings
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : search_strings
    df = filtered_dfview(interactions_dfs[dataset], my_search_strings, min_reads, max_interactions)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df, srna_type, gene_name_types[dataset])
    circos_output = circos_data(df)
    return table_output, cytoscape_output, circos_output, "$(nrow(interactions_dfs[dataset])) $(nrow(df)) $search_strings $my_search_strings"
end