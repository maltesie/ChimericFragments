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
    Input("dropdown-update-layout", "value"),
    Input("ligation", "value"),
    Input("exclusive-search", "value"),
]
const update_selection_states = [
    State("dropdown-update-dataset", "value"),
    State("data-tabs", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, sRNA_type::String, param_dict::Vector{Pair{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, max_fdr, max_bp_fdr, search_strings, layout_value, ligation, exclusive, dataset, tab_value
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    any(isnothing(v) for v in (min_reads, max_interactions, max_bp_fdr, max_fdr)) && throw(PreventUpdate())
    df = filtered_dfview(interactions[dataset], my_search_strings, min_reads, max_interactions, max_fdr, max_bp_fdr, "ligation" in ligation, "exclusive" in exclusive)
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
    Output("plot1", "figure"),
]
update_dataset_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, min_reads::Int) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return min_reads, [Dict("label"=>k, "value"=>k) for k in sort(interactions[dataset].nodes.name)], (
        data = [
            (x = [1, 2, 3], y = [4, 1, 2], type = "bar", name = "SF"),
            (x = [1, 2, 3], y = [2, 4, 5], type = "bar", name = "Montréal"),
        ],
        layout = (title = "Dash Data Visualization",)
    )
end

normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=0, to_max=100) = Int(floor(to_min + value * (to_max-to_min)))
function gene_arrow_and_means(m1::Int, count1::Int, range1::Int, m2::Int, count2::Int, range2::Int, left::Int, right::Int, cds::Int, isnegative::Bool, isrna1::Bool)
    arrows = [
        m1==0 ? html_div(className="colorbar-arrow mean empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow interaction-left" : "colorbar-arrow interaction-right",
                                            title="$m1", style=Dict("left"=>"$(mapvalue(normalize(m1, left, right, false)))%")),
        m2==0 ? html_div(className="colorbar-arrow mode empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow ligation-left" : "colorbar-arrow ligation-right",
                                            title="$m2", style=Dict("left"=>"$(mapvalue(normalize(m2, left, right, false)))%"))
    ]
    p2 = isnegative ? (cds > 0 ? cds : right)-m2+1 : m2-(cds > 0 ? cds : left)+1
    p1 = isnegative ? (cds > 0 ? cds : right)-m1+1 : m1-(cds > 0 ? cds : left)+1
    p2 <= 0 && (p2 -= 1)
    p1 <= 0 && (p1 -= 1)
    s2 = p2>0 ? " +" : " "
    s1 = p1>0 ? " +" : " "
    return html_div(children=[
        #html_p("$m1, $(mapvalue(normalize(m1, left, right, isnegative))), $m2, $(mapvalue(normalize(m2, left, right, isnegative)))"),
        html_p(m2 == 0 ? "no ligation data." : "most frequent ligation point:$s2$p2 ($count2 reads, $range2 nt)", style=Dict("border-top"=>"1px solid #7fafdf", "margin-bottom"=>"12px", "color"=>"DarkSalmon")),
        html_div(className="horizontal deflate", children=[
            html_p(["$left"], className="cb-left"),
            html_div(id=isnegative ? "pointer-left" : "pointer-right", children=arrows),#className="controls-block horizontal", children=arrows),
            html_p(["$right"], className="cb-right")
        ]),
        html_p(m1 == 0 ? "no alignment end data." : "most frequent alignment end:$s1$p1 ($count1 reads, $range1 nt)", style=Dict("margin-top"=>"3px", "border-bottom"=>"1px solid #7fafdf"))
    ])
end

function edge_info(edge_data::Dash.JSON3.Object, interact::Interactions)
    src, dst = parse(Int, edge_data["source"]), parse(Int, edge_data["target"])
    name1, name2 = interact.nodes[src, :name], interact.nodes[dst, :name]
    ref1, ref2 = interact.nodes[src, :ref], interact.nodes[dst, :ref]
    left1, left2 = interact.nodes[src, :left], interact.nodes[dst, :left]
    right1, right2 = interact.nodes[src, :right], interact.nodes[dst, :right]
    return [html_div(id="edge-info", children=[
        html_p("RNA1: $name1 on $ref1 ($(right1-left1+1)nt)"),
        #gene_arrow_and_means(Int(edge_data["modeint1"]), Int(edge_data["modeintcount1"]), Int(edge_data["modeintrange1"]), i1, Int(edge_data["modeligcount1"]), Int(edge_data["modeligrange1"]), l1, r1, c1, strand1=="-", true),
        #html_p(join(", ", edge_data["ligation_points"])), #join(", ", interact.edgestats[(Int(edge_data["src"]), Int(edge_data["dst"]))])),
        #html_br(),
        html_p("RNA2: $name2 on $ref2 ($(right2-left2+1)nt)"),
        #gene_arrow_and_means(Int(edge_data["modeint2"]), Int(edge_data["modeintcount2"]), Int(edge_data["modeintrange2"]), i2, Int(edge_data["modeligcount2"]), Int(edge_data["modeligrange2"]), l2, r2, c2, strand2=="-", false),
        #html_br(),
        #html_p(children=alnstring, style=Dict("white-space" => "pre", "font-family" => "monospace", "max-width"=>"300px", "overflow"=>"scroll", "padding-bottom"=>"13px")),
        html_p("supporting reads: $(edge_data["interactions"]) ($(sum(values(interact.edgestats[src, dst][3]))) with ligation point)"),

    ])]
end

function edge_figure(edge_data::Dash.JSON3.Object, interact::Interactions)
    src, dst = parse(Int, edge_data["source"]), parse(Int, edge_data["target"])
    name1, name2 = interact.nodes[src, :name], interact.nodes[dst, :name]
    ref1, ref2 = interact.nodes[src, :ref], interact.nodes[dst, :ref]
    left1, left2 = interact.nodes[src, :left], interact.nodes[dst, :left]
    right1, right2 = interact.nodes[src, :right], interact.nodes[dst, :right]
    return (
        data = [
            (x = [first(p) for (p,c) in interact.edgestats[(src,dst)][3]], y = [last(p) for (p,c) in interact.edgestats[(src,dst)][3]], type = "scatter", name = "LP"),
        ],
        layout = (title = "$name1 $name2",)
    )
end

ligation_modes_table(ligation_points::Dash.JSON3.Object; ncolumns=5) =
    return isempty(ligation_points) ?
    html_p("none") :
    html_div(
        children=[
            html_p(title="$d", String(k)[1] === '-' ? "$k: $v" : "+$k: $v") for (k,(v,d)) in sort(collect(ligation_points), by=x->parse(Int, String(x[1])))
        ],
        style=Dict("display"=>"grid", "grid-template-columns"=>repeat("1fr ", ncolumns), "max-height"=>"100px", "overflow"=>"scroll")
    )

function node_info(node_data::Dash.JSON3.Object, interactions::Interactions)
    name, ref, strand, left, right, ints = interactions.nodes[parse(Int, node_data["id"]), [:name, :ref, :strand, :left, :right, :nb_ints]]
    return [html_div(id="edge-info", children=[
        html_p("$name"),
        html_p("On $ref ($strand) from $left to $right $(right-left+1)nt)"),
        html_p("$(node_data["nb_partners"]) partner" * (node_data["nb_partners"]>1 ? "s" : "") * " in the current selection."),
        html_br(),
        html_p("Ligation points as RNA1:"),
        #ligation_modes_table(node_data["lig_as_rna1"]),
        html_br(),
        html_p("Ligation points as RNA2:"),
        #ligation_modes_table(node_data["lig_as_rna2"]),
        html_br(),
        html_div(children=[
            html_p("read counts for $name:"),
            html_p("$ints (selection), $ints (total), $ints (single)"),
        ])
    ])]
end

function circos_description(circos_data::Dash.JSON3.Object)
    data = circos_data["data"]
    html_div([
        html_p("RNA1: $(data["name1"]) on $(data["ref1"]) ($(data["strand1"]))"),
        html_p("feature left: $(data["left1"])"),
        html_p("feature right: $(data["right1"])"),
        html_br(),
        html_p("RNA2: $(data["name2"]) on $(data["ref2"]) ($(data["strand2"]))"),
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
    Output("info-output", "children"),
    Output("plotly-graph", "figure"),
]
const update_selected_element_states = [
    State("dropdown-update-dataset", "value")
]
update_selected_element_callback!(app::Dash.DashApp, genome::Dict{String,BioSequences.LongDNA{4}}, interactions::Dict{String, Interactions},
        check_interaction_distances::Tuple{Int,Int}, bp_parameters::NTuple{6,Int}) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs, update_selected_element_states; prevent_initial_call=true) do node_data, edge_data, circos_data, tab_value, dataset
    
    #if tab_value == "circos"
    #    (isnothing(circos_data) || isempty(circos_data)) && return ["Move your mouse over an interaction in the circos plot to display the corresponding partners."]
    #    return [circos_description(circos_data)]
    #elseif tab_value == "table"
    #    return ["The downloadable version of this table contains additional information, e.g. about ligation points."]
    #elseif tab_value == "graph"
        no_node_data = isnothing(node_data) || isempty(node_data)
        no_edge_data = isnothing(edge_data) || isempty(edge_data)
        no_node_data && no_edge_data && 
            return ["test1"], (
                data = [
                    (x = [1, 2, 3], y = [4, 1, 2], type = "bar", name = "SF"),
                    (x = [1, 2, 3], y = [2, 4, 5], type = "bar", name = "Montréal"),
                ],
                layout = (title = "Dash Data Visualization",)
            )
        no_edge_data && return ["test2"], (
            data = [
                (x = [1, 2, 3], y = [4, 1, 2], type = "bar", name = "SF"),
                (x = [1, 2, 3], y = [2, 4, 5], type = "bar", name = "Montréal"),
            ],
            layout = (title = "Dash Data Visualization",)
        )
        src, dst = parse(Int, edge_data[1]["source"]), parse(Int, edge_data[1]["target"])
        return [join([first(p) for (p,c) in interactions[dataset].edgestats[(src,dst)][2]], ", ")], edge_figure(edge_data[1], interactions[dataset])
    #elseif tab_value == "summary"
    #    return ["Change the dataset or selection criteria to update the summary."]
    #end
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
    State("ligation", "value"),
    State("exclusive-search", "value"),
]
const table_column_names = [:name1, :type1, :name2, :type2, :nb_ints, :nb_multi, :in_libs, :pvalue, :fdr, :pred_pvalue, :pred_fdr,
                            :modeint1, :modelig1, :meanlen1, :nms1, :modeint2, :modelig2, :meanlen2, :nms2]
click_table_button_callback!(app::Dash.DashApp, interactions::Dict{String,Interactions}) =
callback!(app, click_table_button_outputs, click_table_button_inputs, click_table_button_states; prevent_initial_call=true) do clicks, dataset, min_reads, max_interactions, max_fdr, max_bp_fdr, search_strings, ligation, exclusive
    if clicks>0
        interact = interactions[dataset]
        my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
        df = DataFrame(filtered_dfview(interact, my_search_strings, min_reads, max_interactions, max_fdr, max_bp_fdr,"ligation" in ligation, "exclusive" in exclusive))
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