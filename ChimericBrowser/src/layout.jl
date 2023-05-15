headline_layout() = html_div(
    id="headline",
    children=[html_h3("ChimericBrowser")]
)

dataset_control_layout(datasets::Vector{String}) = html_div(
    className="control-element",
    children=[
        html_p("Dataset:"),
        dcc_dropdown(
            id="dropdown-update-dataset",
            value=datasets[1],
            clearable=false,
            options=[
                Dict("label"=>p, "value"=>p) for p in datasets
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
            value="clustered",
            clearable=false,
            options=[
                Dict("label"=>l, "value"=>v) for (l,v) in ["clustered"=>"clustered", "grid"=>"grid"]
            ]
        ),
    ]
)

dataset_and_postioning_layout(datasets::Vector{String}) = html_div(
    className="controls-block horizontal",
    children=[dataset_control_layout(datasets), positioning_control_layout()]
)

reads_selection_layout(min_reads::Int, max_fdr::Float64, max_bp_fdr::Float64) = html_div(
    className="controls-block",
    children=[
        html_div(className="horizontal",
            children=[
                html_div(
                    style=Dict("padding-top"=>"2px"),
                    children=[
                        html_p(id="min-reads-text", children=["minimum # of reads:"], style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="min-reads", type="number", value=min_reads, min=min_reads, style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                ),
                html_div(
                    style=Dict("padding-top"=>"2px", "margin-left"=>"15px"),
                    children=[
                        html_p(id="nb-interactions-text", children=["maximum # of interactions:"], style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-interactions", type="number", value=500, style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                )
            ]),
        html_div(
            className="horizontal",
            style=Dict("padding-top"=>"15px"),
            children=[
                html_div(
                    style=Dict("padding-top"=>"2px"),
                    children=[
                        html_p(id="max-fdr-text", children=["max. fisher fdr:"], style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-fdr", type="number", value=max_fdr, min=0.0, max=max_fdr, step=0.01, style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                ),
                html_div(
                    style=Dict("padding-top"=>"2px", "margin-left"=>"45px"),
                    children=[
                        html_p(id="max-bp-fdr-text", children=["max. bp fdr:"], style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-bp-fdr", type="number", value=max_bp_fdr, min=0.0, max=max_bp_fdr, step=0.01, style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                )
            ]),
        dcc_checklist(
            id="ligation",
            style=Dict("padding-top"=>"15px"),
            options = [Dict("label" => "include interactions without ligation data", "value" => "ligation")],
            value = []
        )
    ]
)

search_layout() = html_div(
    className="controls-block",
    children=[
        html_div(
            className="horizontal",
            children=[
                dcc_checklist(
                    id="exclusive-search",
                    options=[Dict("label"=>"exclusive search", "value"=>"exclusive")],
                    style=Dict("width"=>"100%"),
                    value=["exclusive"]
                ),
                dcc_dropdown(
                    id="type-multi-select",
                    options=[],
                    style=Dict("padding-left"=>"10px", "width"=>"100%"),
                    multi=true,
                    placeholder="Select type(s)",
                )
            ]
        ),
        html_div(
            className="horizontal",
            style=Dict("padding-top"=>"10px"),
            children=[
                dcc_dropdown(
                    id="gene-multi-select",
                    options=[],
                    multi=true,
                    placeholder="Select annotation(s)",
                ),
                html_div(style=Dict("width"=>"10px")),
                html_button(
                    "<-",
                    id="add-selected-btn",
                    n_clicks=0,
                    style=Dict("position"=>"relative", "padding"=>"0px", "width"=>"50px", "height"=>"36px", "line-height"=>"0px", "color"=>"#7fafdf"))
            ]
        ),
    ]
)

annotation_control_layout(function_categories::Vector{Dict{String,String}}) = html_div(
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
    id="plots-block",
    children=[
        dcc_graph(id="plotly-graph"),
        html_div(
            className="horizontal",
            children=[
                dcc_clipboard(id="clip"),
                html_p(id="clip-text", "(select a data point)", style=Dict("padding-left"=>"10px")),
            ]
        )

    ]
)

control_column_layout(datasets::Vector{String}, min_reads::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64) = html_div(
    id="left-column",
    children=[
        headline_layout(),
        html_div(
            className="container",
            children=[
                dataset_and_postioning_layout(datasets),
                reads_selection_layout(min_reads, max_fisher_fdr, max_bp_fdr),
                search_layout(),
                info_area_layout()
            ]
        )
    ]
)

cytoscape_layout(stylesheet::Vector{Dict{String,Any}}) = html_div(
    id="graph-container",
    className="container",
    children=[
        cytoscape(
            id="graph",
            elements=[],
            autoRefreshLayout=true,
            stylesheet=stylesheet,
            responsive=true,
            layout=Dict("name"=>"preset", "animate"=>false),
            minZoom=0.2,
            maxZoom=0.9,
            zoom=0.9
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
)

chords_track(data::Vector{Dict{String,Any}}) = Dict(
    "type"=>"CHORDS",
    "data"=>data,
    "config"=>Dict(
        "opacity"=>0.9,
        "color"=>"black",
        "tooltipContent"=>Dict(
        #    "source"=>"RNA1",
        #    "target"=>"RNA2",
        #    "targetEnd"=>"nbints"
        )
    )
)
circos_layout(genome_info::Vector{Pair{String,Int}}) = html_div(
    id="circos-container",
    className="container",
    children=[

        html_div(className="plot", children=[
            dcc_dropdown(id="fdr-source", options=[
                (label = "node degree distribution", value = "degree"),
                (label = "Annotation interaction heatmap", value = "annotation"),
                (label = "Odds ratio distribution", value = "odds"),
                (label = "Basepairing alignments clipping distribution", value = "bp"),
            ], value="bp", clearable=false, multi=false),
            dcc_graph(id="plot1"),
        ]),

        html_div(className="plot", children=[
            dcc_graph(id="plot2"),
        ]),

        html_div(className="plot", children=[
            dcc_slider(id="fdr-value", min=0.0, max=1.0, step=0.01, value=0.1),
            dcc_graph(id="plot3"),
        ]),

        html_div(
            className="plot", children=[
                circos(
                    id="my-dashbio-circos",
                    enableZoomPan=false,
                    enableDownloadSVG=true,
                    size=350,
                    config = Dict(
                        "gap"=>0.003,
                        "cornerRadius"=>5,
                        "innerRadius"=>120,
                        "outerRadius"=>125,
                        "labelRadius"=>140,
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
                    tracks=[chords_track(Dict{String,Any}[])]
                ),
            ]

        ),
    ]
)

table_layout() = html_div(
    id="table-container",
    className="container",
    children=[
        dash_datatable(
            id="table",
            columns=[n == "fdr" ? Dict("name"=>n, "id"=>n, "format"=>Dict("specifier"=>".3f")) : Dict("name"=>n, "id"=>n)
                for n in ["name1", "type1", "name2", "type2", "nb_ints", "fdr", "odds_ratio", "pred_fdr", "in_libs"]],
            style_cell=Dict(
                "height"=> "auto",
                "minWidth"=> "100px", "width"=> "100px", "maxWidth"=> "160px",
                "whiteSpace"=> "normal"
            ),
            style_data_conditional=[
                Dict(
                    "if"=> Dict("state"=> "selected"),
                    "backgroundColor"=> "inherit !important",
                    "border"=> "inherit !important",
                )
            ],

        ),
        html_div([
            html_button("DOWNLOAD CSV", id="btn-csv"),
            dcc_download(id="download-dataframe-csv"),
        ])
    ]
)

summary_layout() = html_div(
    id="summary-container",
    className="container",
    children=[
        html_h3("loading...")
    ]
)

astab(div_component::Component, tab_id::String) = dcc_tab(
    id=lowercasefirst(tab_id) * "-tab",
    className="custom-tab",
    selected_className="custom-tab--selected",
    label=uppercasefirst(tab_id),
    value=lowercasefirst(tab_id),
    children=[div_component]
)

tabs_layout(genome_info::Vector{Pair{String,Int}}, stylesheet::Vector{Dict{String,Any}}) = html_div(
    id="tabs-container",
    className="container",
    children=[
        dcc_tabs(
            id="data-tabs",
            value="summary",
            persistence=true,
            children=[
                astab(cytoscape_layout(stylesheet), "graph"),
                astab(table_layout(), "table"),
                astab(circos_layout(genome_info), "plots"),
                astab(summary_layout(), "summary"),
            ]
        )
    ]
)

browser_layout(datasets::Vector{String}, genome_info::Vector{Pair{String,Int}}, stylesheet::Vector{Dict{String,Any}},
                min_reads::Int, max_fdr::Float64, max_bp_fdr::Float64) = html_div(
    id="root",
    children=[html_div(
        id="app-container",
        children=[
            control_column_layout(datasets, min_reads, max_fdr, max_bp_fdr),
            tabs_layout(genome_info, stylesheet)
        ]
    )]
)