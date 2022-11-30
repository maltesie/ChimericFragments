const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks")
]
const update_selection_inputs = [
    Input("min-reads", "value"),
    Input("max-interactions", "value"),
    Input("gene-multi-select", "value"),
    Input("dropdown-update-layout", "value")
]
const update_selection_states = [
    State("dropdown-update-dataset", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, gene_name_info::Dict{String, Dict{String, Tuple{String, String, Int, Int, Char}}},
    gene_name_position::Dict{String, Dict{String, Dict{String, Float64}}}, sRNA_type::String) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, search_strings, layout_value, dataset
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    df = filtered_dfview(interactions[dataset].edges, my_search_strings, min_reads, max_interactions)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df, gene_name_info[dataset], gene_name_position[dataset], sRNA_type, layout_value)
    circos_output = circos_data(df)
    return table_output, cytoscape_output, circos_output
end

const update_dataset_inputs = [
    Input("dropdown-update-dataset", "value")
]
const update_dataset_outputs = [
    Output("min-reads", "value"),
    Output("gene-multi-select", "options")
]
update_dataset_callback!(app::Dash.DashApp, gene_name_info::Dict{String,Dict{String, Tuple{String, String, Int, Int, Char}}}) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return 1, [Dict("label"=>k, "value"=>k) for k in sort(collect(keys(gene_name_info[dataset])))]
end

normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=0, to_max=100) = Int(floor(to_min + value * (to_max-to_min)))
function css_gradient_and_means(m1::Int, m2::Int, left::Int, right::Int, isnegative::Bool, isrna1::Bool)
    arrows = [
        m1==0 ? html_div(className="colorbar-arrow mean empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow interaction-left" : "colorbar-arrow interaction-right",
                                            style=Dict("left"=>"$(mapvalue(normalize(m1, left, right, false)))%")),
        m2==0 ? html_div(className="colorbar-arrow mode empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow ligation-left" : "colorbar-arrow ligation-right",
                                            style=Dict("left"=>"$(mapvalue(normalize(m2, left, right, false)))%"))
    ]
    return html_div(children=[
        #html_p("$m1, $(mapvalue(normalize(m1, left, right, isnegative))), $m2, $(mapvalue(normalize(m2, left, right, isnegative)))"),
        html_p(m2 == 0 ? "no ligation data." : "ligation points mode: $m2", style=Dict("border-top"=>"1px solid #7fafdf", "margin-bottom"=>"12px", "color"=>"DarkSalmon")),
        html_div(className="horizontal deflate", children=[
            html_p(["$left"], className="cb-left"),
            html_div(id=isnegative ? "pointer-left" : "pointer-right", children=arrows),#className="controls-block horizontal", children=arrows),
            html_p(["$right"], className="cb-right")
        ]),
        html_p(m1 == 0 ? "no interaction data." : "interaction points mode: $m1", style=Dict("margin-top"=>"3px", "border-bottom"=>"1px solid #7fafdf"))
    ])
end

function edge_info(edge_data::Dash.JSON3.Object)
    return [html_div(id="edge-info", children=[
        html_p("RNA1: $(edge_data["source"]) on $(edge_data["ref1"]) ($(edge_data["strand1"]))"),
        css_gradient_and_means(Int(edge_data["modeint1"]), Int(edge_data["modelig1"]), edge_data["left1"], edge_data["right1"], edge_data["strand1"]=="-", true),
        #css_gradient_and_means(edge_data["left1"], edge_data["right1"], edge_data["left1"], edge_data["right1"], edge_data["strand1"]=="-", true),
        html_br(),
        html_p("RNA2: $(edge_data["target"]) on $(edge_data["ref2"]) ($(edge_data["strand2"]))"),
        css_gradient_and_means(Int(edge_data["modeint2"]), Int(edge_data["modelig2"]), edge_data["left2"], edge_data["right2"], edge_data["strand2"]=="-", false),
        #css_gradient_and_means(edge_data["left2"], edge_data["right2"], edge_data["left2"], edge_data["right2"], edge_data["strand2"]=="-", false),
        html_br(),
        html_p("total reads: $(edge_data["interactions"])")
    ])]
end

ligation_modes_table(ligation_points::Dash.JSON3.Object) =
    return isempty(ligation_points) ?
    html_p("none") :
    html_div([
        html_tr([
            html_td("$k:"),
            html_td("$v"),
        ]) for (k,v) in sort(collect(ligation_points), by=x->x[2], rev=true)
    ])
function ligation_modes_arrow_with_tooltip(ligation_points::Dash.JSON3.Object)
end
function node_info(node_data::Dash.JSON3.Object)
    return [html_div(id="edge-info", children=[
        html_p("$(node_data["id"]) on $(node_data["ref"]) ($(node_data["strand"])) has $(node_data["nb_partners"]) partner" * (node_data["nb_partners"]>1 ? "s" : "") * " in the current selection."),
        html_br(),
        html_p("Ligation point modes as RNA1:"),
        ligation_modes_table(node_data["lig_as_rna1"]),
        html_br(),
        html_p("Ligation point modes as RNA2:"),
        ligation_modes_table(node_data["lig_as_rna2"]),
        html_br(),
        html_p("interactions: $(node_data["interactions"])")
    ])]
end

function circos_description(circos_data::Dash.JSON3.Object)
    data = circos_data["data"]
    html_div([
        html_p("RNA1: $(data["source"]) on $(data["ref1"]) ($(data["strand1"]))"),
        html_p("feature left: $(data["left1"])"),
        html_p("feature right: $(data["right1"])"),
        html_br(),
        html_p("RNA2: $(data["target"]) on $(data["ref2"]) ($(data["strand2"]))"),
        html_p("feature left: $(data["left2"])"),
        html_p("feature right: $(data["right2"])"),
        html_br(),
        html_p("interactions: $(data["interactions"])")
    ])
end

const update_selected_element_inputs = [
    Input("graph", "selectedNodeData"),
    Input("graph", "selectedEdgeData"),
    Input("my-dashbio-circos", "eventDatum"),
    Input("data-tabs", "value")
]
const update_selected_element_outputs = [
    Output("info-output", "children")
]
update_selected_element_callback!(app::Dash.DashApp) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs; prevent_initial_call=true) do node_data, edge_data, circos_data, tab_value
    if tab_value == "circos"
        (isnothing(circos_data) || isempty(circos_data)) && return ["Move your mouse over an interaction in the circos plot to display the corresponding partners."]
        return [circos_description(circos_data)]
    elseif tab_value == "table"
        return ["The downloadable version of this table contains additional information, e.g. about ligation points."]
    elseif tab_value == "graph"
        no_node_data = isnothing(node_data) || isempty(node_data)
        no_edge_data = isnothing(edge_data) || isempty(edge_data)
        no_node_data && no_edge_data && return ["Select an edge or node in the graph to display additional information."]
        no_edge_data && return node_info(node_data[1])
        return edge_info(edge_data[1])
    elseif tab_value == "summary"
        return ["Change the dataset or selection criteria to update the summary."]
    end
end

#update_layout_callback!(app::Dash.DashApp) =
#callback!(app, Output("graph", "layout"), Input("dropdown-update-layout", "value"); prevent_initial_call=true) do layout_value
#    return Dict("name"=>layout_value == "random" ? "preset" : layout_value, "animate"=>false)
#end

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
    State("radio-options-csv", "value"),
    State("min-reads", "value"),
    State("max-interactions", "value"),
    State("gene-multi-select", "value"),
]
const table_column_names = [:name1, :type1, :ref1, :strand1, :name2, :type2, :ref2, :strand2, :nb_ints, :nb_multi, :in_libs, :p_value, :fdr,:left1, :right1, :modeint1,
:rel_int1, :modelig1, :rel_lig1, :meanlen1, :nms1, :left2, :right2, :modeint2, :rel_int2, :modelig2, :rel_lig2, :meanlen2, :nms2]
click_table_button_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}) =
callback!(app, click_table_button_outputs, click_table_button_inputs, click_table_button_states; prevent_initial_call=true) do clicks, dataset, csv_option, min_reads, max_interactions, search_strings
    if clicks>0
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        df = (csv_option == "full") ?
            interactions[dataset].edges[!, table_column_names] :
            DataFrame(filtered_dfview(interactions[dataset].edges, my_search_strings, min_reads, max_interactions))[!, table_column_names]
        csvrowwriteriterator = CSV.RowWriter(df)
        dfstring = join(collect(csvrowwriteriterator))
        return [Dict("filename"=>"$(dataset)_$(csv_option)_table.csv", "content"=>dfstring ,"base64"=>false)]
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
        node_data[1]["id"] in my_search_strings && throw(PreventUpdate())
        return push!(my_search_strings, node_data[1]["id"])
    else
        throw(PreventUpdate())
    end
end