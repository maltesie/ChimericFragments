const update_selection_outputs = [
    Output("table", "data"),
    Output("graph", "elements"),
    Output("my-dashbio-circos", "tracks")
]
const update_selection_inputs = [
    Input("min-reads", "value"),
    Input("max-interactions", "value"),
    Input("gene-multi-select", "value"),
]
const update_selection_states = [
    State("dropdown-update-dataset", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, gene_name_type::Dict{String, Dict{String, String}},
    gene_name_position::Dict{String, Dict{String, Dict{String, Float64}}}, sRNA_type::String) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, search_strings, dataset
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    df = filtered_dfview(interactions[dataset].edges, my_search_strings, min_reads, max_interactions)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df, gene_name_type[dataset], gene_name_position[dataset], sRNA_type)
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
update_dataset_callback!(app::Dash.DashApp, gene_name_type::Dict{String,Dict{String,String}}) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return 1, [Dict("label"=>k, "value"=>k) for k in sort(collect(keys(gene_name_type[dataset])))]
end

css_table(data::Vector{Int}, id::String) = html_ul(
    id=id, className="chart", children=[
        html_li(html_span(style=Dict("height"=>"$(Int(floor(v/maximum(data)*100)))%"))) for v in data
    ]
)
normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=-12, to_max=104) = Int(floor(to_min + value * (to_max-to_min)))
function css_gradient_and_means(m1::Int, m2::Int, left::Int, right::Int, isnegative::Bool, id::String)
    arrows = [html_p(id, className="cb-middle")]
    push!(arrows, m1==0 ? html_div(style=Dict("margin-top"=>"-27px")) :
        html_div(className="colorbar-arrow-mean", style=Dict("left"=>"$(mapvalue(normalize(m1, left, right, isnegative)))%")))
    push!(arrows, m2==0 ? html_div(style=Dict("margin-top"=>"-27px")) :
        html_div(className="colorbar-arrow-mode", style=Dict("left"=>"$(mapvalue(normalize(m2, left, right, isnegative)))%")))
    return html_div(children=[
        html_div(className="horizontal deflate", children=[
            html_p(["5'", html_br(), isnegative ? "$right" : "$left"], className="cb-left"),
            html_div(id="colorbar-edges", className="controls-block", children=arrows),
            html_p(["3'", html_br(), isnegative ? "$left" : "$right"], className="cb-right")
        ]),
        m2 == 0 ? html_p("mode of interaction point distribution: $m1") : html_p("mode of ligation point distribution: $m2")
    ])
end
function edge_info(edge_data::Dash.JSON3.Object)
    return [html_div(id="edge-info", children=[
        html_p("$(edge_data["source"])->$(edge_data["target"]): $(edge_data["interactions"]) chimeric read" * (edge_data["interactions"]>1 ? "s" : "") * "."),
        html_br(),
        css_gradient_and_means(Int(edge_data["modeint1"]), Int(edge_data["modelig1"]), edge_data["left1"], edge_data["right1"], edge_data["strand1"]=="-", edge_data["source"]),
        html_br(),
        css_gradient_and_means(Int(edge_data["modeint2"]), Int(edge_data["modelig2"]), edge_data["left2"], edge_data["right2"], edge_data["strand2"]=="-", edge_data["target"]),
    ])]
end

const update_selected_element_inputs = [
    Input("graph", "selectedNodeData"),
    Input("graph", "selectedEdgeData")
]
const update_selected_element_outputs = [
    Output("info-output", "children")
]
#const update_selected_element_states = [
#    State("dropdown-update-dataset", "value")
#]
update_selected_element_callback!(app::Dash.DashApp) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs; prevent_initial_call=true) do node_data, edge_data
    no_node_data = isnothing(node_data) || isempty(node_data)
    no_edge_data = isnothing(edge_data) || isempty(edge_data)
    no_node_data && no_edge_data && return ["Select an edge or node in the graph to display additional information."]
    no_edge_data && return ["$(node_data[1]["id"]) has $(node_data[1]["nb_partners"]) partner" * (node_data[1]["nb_partners"]>1 ? "s" : "") * " in the current selection."]
    return edge_info(edge_data[1])
end