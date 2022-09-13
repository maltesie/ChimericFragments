headline_layout() = html_div(
    id="headline",
    children=[html_h3("ChimericBrowser")]
)

dataset_control_layout(dataset_paths::Vector{String}) = html_div(
    className="control-element",
    children=[
        html_p("Dataset:"),
        dcc_dropdown(
            id="dropdown-update-dataset",
            value=dataset_paths[1],
            clearable=false,
            options=[
                Dict("label"=>basename(p)[1:end-4], "value"=>p) for p in dataset_paths
            ]
        ),
    ]
)

positioning_control_layout() = html_div(
    className="control-element",
    children=[
        html_p("Positioning:"),
        dcc_dropdown(
            id="dropdown-update-layout",
            value="random",
            clearable=false,
            options=[
                Dict("label"=>n, "value"=>n) for n in ["cose", "grid", "random", "circle", "concentric"]
            ]
        ),
    ]
)

reads_selection_layout() = html_div(
    className="controls-block",
    children=[
        html_div(
            className="horizontal",
            children=[
                html_p("Minimum # of reads:"),
                dcc_slider(
                    id="reads-slider",
                    min=0,
                    max=100,
                    step=1,
                    value=50
                ),
            ]
        ),
        html_div(
            className="horizontal",
            style=Dict("padding-top"=>"10px"),
            children=[
                html_p(id="nb-interactions-text", style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                dcc_input(id="max-interactions", type="number", value=300, style=Dict("max-height"=>"18px", "min-width"=>"80px"))
            ]
        )
    ]
)

search_layout() = html_div(
    className="controls-block",
    children=[
        html_p("Search:"),
        html_div(
            className="horizontal",
            children=[
                dcc_dropdown(
                    id="gene-multi-select",
                    options=[],
                    multi=true
                ),
                html_div(style={"width"=>"10px"}),
                html_button(
                    "<-",
                    id="add-selected-btn",
                    n_clicks=0,
                    style=Dict("position"=>"relative", "padding"=>"0px", "width"=>"50px", "height"=>"36px", "line-height"=>"0px", "color"=>"#7fafdf"))
            ]
        )
    ]
)

annotation_control_layout(function_categories::Voctor{Dict{String,String}}) = html_div(
    className="controls-block",
    children=[
        html_p("Annotated Functions:"),
        dcc_dropdown(
            id="function-multi-select",
            options=function_categories,
            clearable=false,
            multi=true
        ),
        dcc_checklist(
            id="color-checklist",
            options=[Dict("label"=>"Restrict current graph to selected functions", "value"=>"restrict_fun")],
            style=Dict("padding-top"=>"12px", "display"=>"flex", "flex-direction"=>"column")
        )
    ]
)

info_area_layout() = html_div(
    className="controls-block",
    children=[html_p(id="info-output", children=["No interaction selected."])]
)

cytoscape_layout() = html_div(
    id="graph-container",
    className="container",
    children=[
        cytoscape(
            id="graph",
            elements=[],
            autoRefreshLayout=True,
            stylesheet=stylesheet,
            responsive=True,
            layout=Dict("name"=>"random"),
            minZoom=0.1,
            maxZoom=2.0,
        ),
        html_div(
            children=[
                html_button("DOWNLOAD SVG", id="save-svg", n_clicks=0),
                html_div(style=Dict("height"=>"8%")),
                html_div(id="legend-container"),
                html_div(style=Dict("height"=>"100%"))
            ]
        )
    ]
),

circos_layout(genome_info::Vector{Pair{String,Int}}) = html_div(
    id="circos-container",
    className="container",
    children=[
        circos(
            id="circos-simple",
            config = Dict(
                "gap"=>0.003,
                "cornerRadius"=>5,
                "innerRadius"=>290,
                "ticks"=> Dict(
                    "display"=>true,
                    "spacing"=>100000,
                    "color"=>"#000",
                    "labelDenominator"=>1000000,
                    "labelSuffix"=>" Mb",
                    "labelFont"=>"Arial",
                    "majorSpacing"=>5,
                    "minorSpacing"=>1,
                    "labelSpacing"=>5
                )
            ),
            layout=[Dict("id"=> n, "label"=> n, "len"=> l) for (n,l) in genome_info],
            selectEvent=Dict("2"=> "both"),
            tracks=[Dict(
                "type"=>"CHORDS",
                "data"=>[],
                "config"=>Dict(
                    "opacity"=>0.9,
                    "color"=>Dict("name"=> "color"),
                    "tooltipContent"=>Dict("name":"nb_ints")
                )
            )]
        )
    ]
)

table_layout() = html_div(
    id="table-container",
    className="container",
    children=[
        dash_datatable(
            id="table",
            columns=[Dict("name"=> n, "id"=> n) for n in ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]],
            style_cell=Dict(
                "height"=> "auto",
                "minWidth"=> "140px", "width"=> "140px", "maxWidth"=> "140px",
                "whiteSpace"=> "normal"
            ),
            style_data_conditional=[
                Dict(
                    "if"=> Dict("state"=> "selected"),
                    "backgroundColor"=> "inherit !important",
                    "border"=> "inherit !important",
                )
            ]
        ),
        html_div([
            html_button("DOWNLOAD CSV", id="btn_csv"),
            dcc_download(id="download-dataframe-csv"),
        ])
    ]
)

browser_layout = html_div([circos_layout, cyto_layout])