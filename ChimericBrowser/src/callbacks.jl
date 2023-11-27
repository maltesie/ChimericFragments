# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks"),
    Output("summary-container", "children"),
    Output("data-tabs", "value"),
    Output("graph", "selectedNodeData"),
    Output("graph", "selectedEdgeData"),
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
    State("graph", "selectedEdgeData"),
]
# Callbacks get executed when input variables change their value. update_selection_callback runs when any of the search or filter criteria
# in the control elements change. It first filters the total amount of data in the selected dataset and computes inputs for the different
# visualization modes (graph with cytoscape, circos plot, table view and dataset summary)
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, param_dict::Vector{Pair{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, max_fdr,
        max_bp_fdr, search_strings, type_strings, layout_value, ligation, exclusive, dataset, tab_value, selected_node, selected_edge

    # generate list of genes to search for
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)

    # generate list of annotation types to search for
    my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)

    # handle misformatted inputs for float values
    any(isnothing(v) for v in (min_reads, max_interactions, max_bp_fdr, max_fdr)) && throw(PreventUpdate())

    # create DataFrame view according to filter and search criteria
    df = filtered_dfview(interactions[dataset], my_search_strings, my_type_strings, min_reads, max_interactions, Float64(max_fdr),
        Float64(max_bp_fdr), "ligation" in ligation, "exclusive" in exclusive)

    # compute inputs for visualization modes (graph with cytoscape, circos plot, table view and dataset summary)
    table_output = table_data(df, interactions[dataset])
    cytoscape_output = cytoscape_elements(df, interactions[dataset], layout_value)
    circos_output = circos_data(df, interactions[dataset])
    summary_output = summary_statistics(df, interactions[dataset], param_dict)
    return table_output, cytoscape_output, circos_output, summary_output, tab_value, selected_node, selected_edge
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
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
# Callbacks get executed when input variables change their value. update_dataset_callback runs when another dataset is selected
update_dataset_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, min_reads::Int) =
callback!(app, update_dataset_outputs, update_dataset_inputs, update_dataset_states; prevent_initial_call=false) do dataset, fdr
    # Collect genes and annotation types present in the currently selected dataset
    return min_reads,
    [Dict("label"=>k, "value"=>k) for k in sort(interactions[dataset].nodes.name)],
    [Dict("label"=>k, "value"=>k) for k in unique(interactions[dataset].nodes.type)],
    fdr
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
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
# Callbacks get executed when input variables change their value. update_plots_callback runs when a new plot type is selected or
# the fdr-value slider changes. The circos pot in the plots area is controlled by update_selection_callback.
update_plots_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, randseq_model_ecdf::ECDF, genome_model_ecdf::ECDF) =
callback!(app, update_plots_outputs, update_plots_inputs, update_plots_states; prevent_initial_call=false) do plot_type, plot_fdr, dataset
    # generate plots for the selected plot type
    p1, p2 = plot_pair(interactions[dataset], plot_type, Float64(plot_fdr))
    # the score distribution plot is always displayed and chages according to the selected FDR value.
    return bp_score_dist_plot(interactions[dataset], randseq_model_ecdf, genome_model_ecdf, Float64(plot_fdr)), p1, p2
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
const update_aggregation_slider_inputs = [
    Input("aggregation-fdr-slider", "value")
]
const update_aggregation_slider_outputs = [
    Output("plotly-graph", "figure")
]
const update_aggregation_slider_states = [
    State("graph", "selectedNodeData"),
    State("graph", "selectedEdgeData"),
    State("dropdown-update-dataset", "value"),
]
# Callbacks get executed when input variables change their value. update_aggregation_slider_callback runs when the slider for the local FDR
# changes and affects the ligation points and aggregated basepairing plots.
update_aggregation_slider_callback!(app::Dash.DashApp, genome::Dict{String,BioSequences.LongDNA{4}}, interactions::Dict{String, Interactions},
        check_interaction_distances::Tuple{Int,Int}, model::AffineGapScoreModel) =
callback!(app, update_aggregation_slider_outputs, update_aggregation_slider_inputs, update_aggregation_slider_states; prevent_initial_call=true) do aggregation_fdr, node_data, edge_data, dataset

    # check if a node or an edge is currently selected
    no_node_data = isnothing(node_data) || isempty(node_data)
    no_edge_data = isnothing(edge_data) || isempty(edge_data)

    # return static empty plot with hint to select node or edge as title
    no_node_data && no_edge_data && return [empty_figure]

    # return baspair aggregation plot if a node is selected. node_figure is defined in plots.jl
    no_edge_data && return [node_figure(node_data[1], interactions[dataset], check_interaction_distances, Float64(aggregation_fdr))]

    # if an edge is selected, return ligation points plot. edge_figure is defined in plots.jl
    return [edge_figure(edge_data[1], interactions[dataset], genome, check_interaction_distances, model, Float64(aggregation_fdr))]
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
const update_selected_element_inputs = [
    Input("graph", "selectedNodeData"),
    Input("graph", "selectedEdgeData"),
]
const update_selected_element_outputs = [
    Output("aggregation-fdr-slider", "value"),
]
const update_selected_element_states = [
    State("aggregation-fdr-slider", "value"),
]
# Callbacks get executed when input variables change their value. update_selected_element_callback is an axuiliary callback which sets the
# local fdr to its current value again to trigger new baspairing plots when the selected node or edge changes.
update_selected_element_callback!(app::Dash.DashApp) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs, update_selected_element_states; prevent_initial_call=true) do _, _, aggregation_fdr
    return [aggregation_fdr]
end

# Callback for downloading the current graph as a .svg
click_cyto_button_callback!(app::Dash.DashApp) =
callback!(app, Output("graph", "generateImage"), Input("save-svg", "n_clicks"), State("dropdown-update-dataset", "value"); prevent_initial_call=true) do clicks, dataset
    clicks>0 && return Dict("type"=>"svg", "action"=>"download", "filename"=>"$(dataset)_graph")
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
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
# sets the order in which coulumns in the downloaded table appear.
const table_column_names = [:name1, :type1, :name2, :type2, :nb_ints, :nb_multi, :in_libs,
                            :fisher_pvalue, :fisher_fdr, :odds_ratio, :bp_pvalue, :bp_fdr, :meanlen1, :nms1, :meanlen2, :nms2]
# Callback for downloading the currently selected interactions as a .csv file
click_table_button_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}) =
callback!(app, click_table_button_outputs, click_table_button_inputs, click_table_button_states; prevent_initial_call=true) do clicks, dataset,
        min_reads, max_interactions, max_fisher_fdr, max_bp_fdr, search_strings, type_strings, ligation, exclusive
    if clicks>0
        # perform same filtering as in update_selection_callback, first select current dataset
        interact = interactions[dataset]
        # generate list of genes to search for
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        # generate list of annotation types to search for
        my_type_strings = isnothing(type_strings) || all(isempty.(type_strings)) ? String[] : string.(type_strings)
        # create DataFrame view according to filter and search criteria
        df = DataFrame(filtered_dfview(interact, my_search_strings, my_type_strings, min_reads, max_interactions,
            Float64(max_fisher_fdr), Float64(max_bp_fdr), "ligation" in ligation, "exclusive" in exclusive))
        # collect gene names and annotation types
        df.name1 = interact.nodes.name[df.src]
        df.name2 = interact.nodes.name[df.dst]
        df.type1 = interact.nodes.type[df.src]
        df.type2 = interact.nodes.type[df.dst]
        # select coulms in above defined order
        df = df[!, table_column_names]
        # create csv writer and write file into temporary string
        csvrowwriteriterator = CSV.RowWriter(df)
        dfstring = join(collect(csvrowwriteriterator))
        # use downloader component from Dash to send the .csv to the browser
        return [Dict("filename"=>"$(dataset)_table.csv", "content"=>dfstring ,"base64"=>false)]
    end
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
const click_add_node_inputs = [
    Input("add-selected-btn", "n_clicks")
]
const click_add_node_outputs = Output("gene-multi-select", "value")

const click_add_node_states = [
    State("graph", "selectedNodeData"),
    State("gene-multi-select", "value"),
]
# Callbacks get executed when input variables change their value. update_add_node_callback handles the arrow button that adds the
# currently selected node to the list of genes to search for.
click_add_node_callback!(app::Dash.DashApp) =
callback!(app, click_add_node_outputs, click_add_node_inputs, click_add_node_states; prevent_initial_call=true) do clicks, node_data, search_strings
    # check, if a node is selected, otherwise do nothing
    if clicks>0 && !isnothing(node_data) && !isempty(node_data)
        # collect currently selected search terms and add currently selected node to them
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        node_data[1]["name"] in my_search_strings && throw(PreventUpdate())
        return push!(my_search_strings, node_data[1]["name"])
    else
        throw(PreventUpdate())
    end
end

# Define input, output and state variables refering to values of components in the graphical interface defined in layout.jl
const click_clipboard_inputs = [
    Input("plotly-graph", "clickData")
]
const click_clipboard_outputs = [
    Output("clip-text", "children"),
    Output("clip", "data"),
]
const click_clipboard_states = [
    State("graph", "selectedNodeData"),
    State("graph", "selectedEdgeData"),
    State("dropdown-update-dataset", "value"),
    State("aggregation-fdr-slider", "value")
]
# Callbacks get executed when input variables change their value. clock_clipboard_callback is executed when a data point in a ligation points plot
# or a basepair aggregation plot is clicked and sends the tooltip text to the clipboard component.
click_clipboard_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}, bp_len::Int) =
callback!(app, click_clipboard_outputs, click_clipboard_inputs, click_clipboard_states; prevent_initial_call=true) do click_data, node_data, edge_data, dataset, max_fdr
    # check if any datapoint is currently selected
    if !isnothing(click_data) && !isempty(click_data)
        # check if the current plot is a basepairing aggregation plot or a ligation point plot
        if !isnothing(node_data) && !isempty(node_data)
            # coordinate of clicked data point
            p = click_data["points"][1]["x"]

            #index of currently selected gene
            idx = parse(Int, node_data[1]["id"])
            n = interactions[dataset].nodes.name[idx]

            # compute counts for all partners for the orientation given by the clicked data point
            counts = click_data["points"][1]["curveNumber"] == 0 ?
                count_ligation_sites_as1(idx, node_data[1], interactions[dataset], bp_len, Float64(max_fdr)) :
                count_ligation_sites_as2(idx, node_data[1], interactions[dataset], bp_len, Float64(max_fdr))

            p in keys(counts) || throw(PreventUpdate())
            # select only partners at clicked position
            names_and_counts = sort(collect(counts[p]), by=x->x[2], rev=true)

            # create DataFrame and csv writer and write content in CSV format into temporary string
            df = DataFrame("name"=>interactions[dataset].nodes.name[first.(names_and_counts)], "count"=>last.(names_and_counts))
            csvrowwriteriterator = CSV.RowWriter(df)
            dfstring = join(collect(csvrowwriteriterator))

            rna_type = click_data["points"][1]["curveNumber"] == 0 ? "RNA1" : "RNA2"

            # return description and initiate download
            "partners of $n at position $p", Dict("filename"=>"$(dataset)_$(n)_$(rna_type)_$(p)_partners.csv", "content"=>dfstring ,"base64"=>false)

        elseif !isnothing(edge_data) && !isempty(edge_data)
            # get indices and names of the interaction
            src, dst = parse(Int, edge_data[1]["source"]), parse(Int, edge_data[1]["target"])
            nsrc, ndst = interactions[dataset].nodes.name[src], interactions[dataset].nodes.name[dst]

            # return description and initiate download of tooltip text containing info and a basepairing prediction for the currently selected ligation point.
            "basepairing between $nsrc and $ndst", Dict("filename"=>"$(dataset)_$(nsrc)_$(ndst)_basepairing.txt", "content"=>"$(replace(click_data["points"][1]["text"], "<br>"=>"\n"))" ,"base64"=>false)
        end
    else
        throw(PreventUpdate())
    end
end