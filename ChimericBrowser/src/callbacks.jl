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
    Input("ligation", "value"),
    Input("gene-multi-select", "value"),
    Input("dropdown-update-layout", "value")
]
const update_selection_states = [
    State("dropdown-update-dataset", "value"),
    State("data-tabs", "value")
]
update_selection_callback!(app::Dash.DashApp, interactions::Dict{String, Interactions}, gene_name_info::Dict{String, Dict{String, Tuple{String, String, Int, Int, Char, Int, Int}}},
    gene_name_position::Dict{String, Dict{String, Dict{String, Float64}}}, sRNA_type::String, param_dict::Vector{Pair{String, String}}) =
callback!(app, update_selection_outputs, update_selection_inputs, update_selection_states; prevent_initial_call=true) do min_reads, max_interactions, max_fdr, max_bp_fdr, ligation, search_strings, layout_value, dataset, tab_value
    my_search_strings = isnothing(search_strings) || all(isempty.(search_strings)) ? String[] : string.(search_strings)
    df = filtered_dfview(interactions[dataset].edges, my_search_strings, min_reads, max_interactions, max_fdr, max_bp_fdr, "ligation" in ligation)
    table_output = table_data(df)
    cytoscape_output = cytoscape_elements(df, interactions[dataset], gene_name_info[dataset], gene_name_position[dataset], sRNA_type, layout_value)
    circos_output = circos_data(df)
    summary_output = summary_statistics(interactions[dataset], df, param_dict)
    return table_output, cytoscape_output, circos_output, summary_output, tab_value
end

const update_dataset_inputs = [
    Input("dropdown-update-dataset", "value")
]
const update_dataset_outputs = [
    Output("min-reads", "value"),
    Output("gene-multi-select", "options")
]
update_dataset_callback!(app::Dash.DashApp, gene_name_info::Dict{String,Dict{String, Tuple{String, String, Int, Int, Char, Int, Int}}}, min_reads::Int) =
callback!(app, update_dataset_outputs, update_dataset_inputs; prevent_initial_call=false) do dataset
    return min_reads, [Dict("label"=>k, "value"=>k) for k in sort(collect(keys(gene_name_info[dataset])))]
end

normalize(value::Int, mi::Int, ma::Int, rev::Bool) = rev ? 1-(value-mi)/(ma-mi) : (value-mi)/(ma-mi)
mapvalue(value::Float64; to_min=0, to_max=100) = Int(floor(to_min + value * (to_max-to_min)))
function gene_arrow_and_means(m1::Int, count1::Int, m2::Int, count2::Int, left::Int, right::Int, isnegative::Bool, isrna1::Bool)
    arrows = [
        m1==0 ? html_div(className="colorbar-arrow mean empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow interaction-left" : "colorbar-arrow interaction-right",
                                            title="$count1 reads", style=Dict("left"=>"$(mapvalue(normalize(m1, left, right, false)))%")),
        m2==0 ? html_div(className="colorbar-arrow mode empty") : html_div(className=isrna1 != isnegative ? "colorbar-arrow ligation-left" : "colorbar-arrow ligation-right",
                                            title="$count2 reads", style=Dict("left"=>"$(mapvalue(normalize(m2, left, right, false)))%"))
    ]
    p2 = isnegative ? right-m2+1 : m2-left+1
    p1 = isnegative ? right-m1+1 : m1-left+1
    s2 = p2>=0 ? " +" : " "
    s1 = p1>=0 ? " +" : " "
    return html_div(children=[
        #html_p("$m1, $(mapvalue(normalize(m1, left, right, isnegative))), $m2, $(mapvalue(normalize(m2, left, right, isnegative)))"),
        html_p(m2 == 0 ? "no ligation data." : "most frequent ligation point:$s2$p2 ($count2 reads)", style=Dict("border-top"=>"1px solid #7fafdf", "margin-bottom"=>"12px", "color"=>"DarkSalmon")),
        html_div(className="horizontal deflate", children=[
            html_p(["$left"], className="cb-left"),
            html_div(id=isnegative ? "pointer-left" : "pointer-right", children=arrows),#className="controls-block horizontal", children=arrows),
            html_p(["$right"], className="cb-right")
        ]),
        html_p(m1 == 0 ? "no alignment end data." : "most frequent alignment end:$s1$p1 ($count1 reads)", style=Dict("margin-top"=>"3px", "border-bottom"=>"1px solid #7fafdf"))
    ])
end

alnchar(x::DNA, y::DNA) =
    if (x == DNA_A && y == DNA_T) || (x == DNA_T && y == DNA_A) || (x == DNA_C && y == DNA_G) || (x == DNA_G && y == DNA_C)
        '|'
    elseif (x == DNA_G && y == DNA_T) || (x == DNA_T && y == DNA_G)
        '⋅'
    else
        ' '
    end
function baisepairing_string(aln::PairwiseAlignment, offset1::Int, offset2::Int; width::Integer=20)
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    # width of position numbers
    posw = ndigits(max(offset1 + anchors[end].seqpos, offset2 - anchors[1].refpos)) + 1
    outstring = ""
    i = 0
    seqpos = offset1 + anchors[1].seqpos
    refpos = offset2 - anchors[1].refpos + 2
    seqbuf = IOBuffer()
    refbuf = IOBuffer()
    matbuf = IOBuffer()
    next_xy = iterate(aln)
    while next_xy !== nothing
        (x, y), s = next_xy
        next_xy = iterate(aln ,s)

        i += 1
        if x != gap(eltype(seq))
            seqpos += 1
        end
        if y != gap(eltype(ref))
            refpos -= 1
        end

        if i % width == 1
            print(seqbuf, " RNA1:", lpad(seqpos, posw), ' ')
            print(refbuf, " RNA2:", lpad(refpos, posw), ' ')
            print(matbuf, " "^(posw + 7))
        end

        print(seqbuf, RNA(x))
        print(refbuf, RNA(y))
        print(matbuf, alnchar(x, y))

        if i % width == 0
            print(seqbuf, lpad(seqpos, posw))
            print(refbuf, lpad(refpos, posw))
            print(matbuf)

            outstring *= String(take!(seqbuf)) * "\n" * String(take!(matbuf)) * "\n" * String(take!(refbuf)) * "\n\n"

            if next_xy !== nothing
                seek(seqbuf, 0)
                seek(matbuf, 0)
                seek(refbuf, 0)
            end
        end
    end

    if i % width != 0
        print(seqbuf, lpad(seqpos, posw))
        print(refbuf, lpad(refpos, posw))
        print(matbuf)

        outstring *= String(take!(seqbuf)) * "\n" * String(take!(matbuf)) * "\n" * String(take!(refbuf)) * "\n\n"
    end
    outstring
end

const scores = Dict((DNA_A, DNA_T)=>4, (DNA_T, DNA_A)=>4, (DNA_C, DNA_G)=>5, (DNA_G, DNA_C)=>5, (DNA_G, DNA_T)=>1, (DNA_T, DNA_G)=>1)
const model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-5, default_mismatch=-5); gap_open=-6, gap_extend=-5)
function edge_info(edge_data::Dash.JSON3.Object, genome::Dict{String, BioSequences.LongDNA{4}}, check_interaction_distances::Tuple{Int,Int})
    i1, i2 = Int(edge_data["modelig1"]), Int(edge_data["modelig2"])
    ref1, ref2 = edge_data["ref1"], edge_data["ref2"]
    strand1, strand2 = edge_data["strand1"], edge_data["strand2"]
    l1, l2, r1, r2 = Int(edge_data["left1"]), Int(edge_data["left2"]), Int(edge_data["right1"]), Int(edge_data["right2"])
    alnstring = if (Int(edge_data["modelig1"]) != 0) && (Int(edge_data["modelig2"]) != 0)
        s1 = strand1=="+" ?
            genome[ref1][(i1-check_interaction_distances[1]):(i1+check_interaction_distances[2])] :
            BioSequences.reverse_complement(genome[ref1][(i1-check_interaction_distances[2]):(i1+check_interaction_distances[1])])
        s2 = strand2=="-" ?
            BioSequences.complement(genome[ref2][(i2-check_interaction_distances[1]):(i2+check_interaction_distances[2])]) :
            BioSequences.reverse(genome[ref2][(i2-check_interaction_distances[2]):(i2+check_interaction_distances[1])])
        p = pairalign(LocalAlignment(), s1, s2, model)
        baisepairing_string(alignment(p),
            (strand1=="+" ? ((i1-check_interaction_distances[1])-l1) : (r1-(i1+check_interaction_distances[1]))),
            (strand2=="-" ? (r2-(i2-check_interaction_distances[1])) : ((i2+check_interaction_distances[1])-l2)))
    else
        "no ligation data."
    end
    return [html_div(id="edge-info", children=[
        html_p("RNA1: $(edge_data["source"]) on $ref1"),
        gene_arrow_and_means(Int(edge_data["modeint1"]), Int(edge_data["modeintcount1"]), i1, Int(edge_data["modeligcount1"]), l1, r1, strand1=="-", true),
        html_br(),
        html_p("RNA2: $(edge_data["target"]) on $ref2"),
        gene_arrow_and_means(Int(edge_data["modeint2"]), Int(edge_data["modeintcount2"]), i2, Int(edge_data["modeligcount2"]), l2, r2, strand2=="-", false),
        html_br(),
        html_p(children=alnstring, style=Dict("white-space" => "pre-wrap", "font-family" => "monospace")),
        html_br(),
        html_p("total reads: $(edge_data["interactions"]) ($(edge_data["ligcount"]) ligation points)")
    ])]
end

ligation_modes_table(ligation_points::Dash.JSON3.Object; ncolumns=5) =
    return isempty(ligation_points) ?
    html_p("none") :
    html_div(
        children=[
            html_p(title="$d", String(k)[1] === '-' ? "$k: $v" : "+$k: $v") for (k,(v,d)) in sort(collect(ligation_points), by=x->parse(Int, String(x[1])))
        ],
        style=Dict("display"=>"grid", "grid-template-columns"=>repeat("1fr ", ncolumns))
    )

function node_info(node_data::Dash.JSON3.Object)
    return [html_div(id="edge-info", children=[
        html_p("$(node_data["id"])"),
        html_p("On $(node_data["ref"]) ($(node_data["strand"])) from $(node_data["left"]) to $(node_data["right"]) ($(node_data["right"]-node_data["left"]+1)nt)"),
        html_p("$(node_data["nb_partners"]) partner" * (node_data["nb_partners"]>1 ? "s" : "") * " in the current selection."),
        html_br(),
        html_p("Ligation points as RNA1:"),
        ligation_modes_table(node_data["lig_as_rna1"]),
        html_br(),
        html_p("Ligation points as RNA2:"),
        ligation_modes_table(node_data["lig_as_rna2"]),
        html_br(),
        html_div(children=[
            html_p("read counts for $(node_data["id"]):"),
            html_p("$(node_data["interactions"]) (selection), $(node_data["nb_ints_total"]) (total), $(node_data["nb_single"]) (single)"),
        ])
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
update_selected_element_callback!(app::Dash.DashApp, genome::Dict{String,BioSequences.LongDNA{4}}, check_interaction_distances::Tuple{Int,Int}) =
callback!(app, update_selected_element_outputs, update_selected_element_inputs; prevent_initial_call=true) do node_data, edge_data, circos_data, tab_value
    if tab_value == "circos"
        (isnothing(circos_data) || isempty(circos_data)) && return ["Move your mouse over an interaction in the circos plot to display the corresponding partners."]
        return [circos_description(circos_data)]
    elseif tab_value == "table"
        return ["The downloadable version of this table contains additional information, e.g. about ligation points."]
    elseif tab_value == "graph"
        no_node_data = isnothing(node_data) || isempty(node_data)
        no_edge_data = isnothing(edge_data) || isempty(edge_data)
        no_node_data && no_edge_data && return ["Select an edge or node in the graph to display additional information. Click the <- button to add a selected node to the search."]
        no_edge_data && return node_info(node_data[1])
        return edge_info(edge_data[1], genome, check_interaction_distances)
    elseif tab_value == "summary"
        return ["Change the dataset or selection criteria to update the summary."]
    end
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