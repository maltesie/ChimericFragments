test_callback!(app::Dash.DashApp) = 
callback!(app, Output("info-output", "children"), Input("reads-slider", "value")) do slider_value
    "$slider_value"
end

const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks")
]

const update_selection_inputs = [
    Input("dropdown-update-dataset", "value"),
    Input("reads-slider", "value"),
    Input("max-interactions", "value"),
    Input("gene-multi-select", "value"),
]

update_selection_callback!(app::Dash.DashApp, interactions_dfs::Dict{String,DataFrame}) =
callback!(app, update_selection_outputs, update_selection_inputs) do dataset, min_reads, max_interactions, search_strings
    df = filtered_view(interactions_dfs[dataset], search_strings, min_reads, max_interactions)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df)
    circos_output = circos_data(df)
    return table_output, cytoscape_output, circos_output 
end