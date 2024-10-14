# layout.jl defines the html structure of the graphical interface. Documentation available at https://dash.plotly.com/julia/layout
# the style attribute of every html element directly sets the CSS attributes of the component.

# Set the main title on the upper left corner
headline_layout() = html_div(
    id="headline",
    children=[html_h3("ChimericBrowser")]
)

# define the dropbox for choosing the current dataset from all available datasets in the jld folder created by analyze.jl
dataset_control_layout(datasets::Vector{String}) = html_div(
    className="control-element", # control-element is defined in assets/general.css
    children=[
        html_p("Dataset:"),
        dcc_dropdown(
            id="dropdown-update-dataset", # used to access the element in the callbacks
            value=datasets[1],
            clearable=false,
            options=[
                Dict("label"=>p, "value"=>p) for p in datasets
            ]
        ),
    ]
)

# define the dropbox for choosing the node positioning algorithm
positioning_control_layout() = html_div(
    className="control-element", # control-element is defined in assets/general.css
    children=[
        html_p("Positioning:"),
        dcc_dropdown(
            id="dropdown-update-layout", # used to access the element in the callbacks
            value="clustered",
            clearable=false,
            options=[
                Dict("label"=>l, "value"=>v) for (l,v) in ["clustered"=>"clustered", "grid"=>"grid"]
            ]
        ),
    ]
)

# combine dataset selection and node positioning next to each other
dataset_and_postioning_layout(datasets::Vector{String}) = html_div(
    className="controls-block horizontal", # controls-block and horizontal are defined in assets/general.css
    children=[dataset_control_layout(datasets), positioning_control_layout()]
)

reads_selection_layout(min_reads::Int, max_fdr::Float64, max_bp_fdr::Float64) = html_div(
    className="controls-block", # controls-block is defined in assets/general.css
    children=[
        html_div(className="horizontal", # horizontal is defined in assets/general.css
            children=[
                # control element for selecting the minimum number of reads per interaction
                html_div(
                    style=Dict("padding-top"=>"2px"),
                    children=[
                        html_p(id="min-reads-text", children=["minimum # of reads:"],
                            style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="min-reads", type="number", value=min_reads, min=min_reads,
                            style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                ),
                # control element for setting the maximum number of interactions in the selection
                html_div(
                    style=Dict("padding-top"=>"2px", "margin-left"=>"15px"),
                    children=[
                        html_p(id="nb-interactions-text", children=["maximum # of interactions:"],
                            style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-interactions", type="number", value=500,
                            style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                )
            ]),
        html_div(
            className="horizontal", # horizontal is defined in assets/general.css
            style=Dict("padding-top"=>"15px"),
            children=[
                # control element for setting the maximum FDR according to Fishers exact test
                html_div(
                    style=Dict("padding-top"=>"2px"),
                    children=[
                        html_p(id="max-fdr-text", children=["max. fisher fdr:"],
                            style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-fdr", type="number", value=max_fdr, min=0.0, max=max_fdr, step=0.01,
                            style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                ),
                # control element for setting the maximum FDR according to complementarity around ligation points
                html_div(
                    style=Dict("padding-top"=>"2px", "margin-left"=>"45px"),
                    children=[
                        html_p(id="max-bp-fdr-text", children=["max. bp fdr:"],
                            style=Dict("padding-right"=>"5px", "padding-top"=>"2px")),
                        dcc_input(id="max-bp-fdr", type="number", value=max_bp_fdr, min=0.0, max=max_bp_fdr, step=0.01,
                            style=Dict("max-height"=>"18px", "min-width"=>"70px", "max-width"=>"70px"))
                    ]
                )
            ]),
        # control element for including interactions without ligation points
        dcc_checklist(
            id="ligation",
            style=Dict("padding-top"=>"15px"),
            options = [Dict("label" => "include interactions without ligation data", "value" => "ligation")],
            value = ["ligation"]
        )
    ]
)

# combine filtering criteria for number of interactions, number of reads per interaction and FDRs
search_layout() = html_div(
    className="controls-block", # controls-block is defined in assets/general.css
    children=[
        html_div(
            className="horizontal", # horizontal is defined in assets/general.css
            children=[
                # control element for toggling AND and OR behavior of combining search terms (genes and annotation types)
                dcc_checklist(
                    id="exclusive-search",
                    options=[Dict("label"=>"exclusive search", "value"=>"exclusive")],
                    style=Dict("width"=>"100%"),
                    value=["exclusive"]
                ),
                # control element for selecting search terms for annotation types
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
            className="horizontal", # horizontal is defined in assets/general.css
            style=Dict("padding-top"=>"10px"),
            children=[
                # control element for selecting search terms for gene names
                dcc_dropdown(
                    id="gene-multi-select",
                    options=[],
                    multi=true,
                    placeholder="Select annotation(s)",
                ),
                html_div(style=Dict("width"=>"10px")),
                # button to add currently selected node to the list of genes to filter for
                html_button(
                    "<-",
                    id="add-selected-btn",
                    n_clicks=0,
                    style=Dict("position"=>"relative", "padding"=>"0px", "width"=>"50px", "height"=>"36px", "line-height"=>"0px", "color"=>"#7fafdf"))
            ]
        ),
    ]
)

# control element for selecting functional annotations, currently not used
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

# html element for displaying ligation points plots and complementarity plots
info_area_layout() = html_div(
    id="plots-block", # plots-block is defined in assets/general.css
    children=[
        dcc_graph(id="plotly-graph"),
        html_div(style=Dict("background-color"=>"white", "display"=>"flex", "flex-direction"=>"row", "padding-left"=>"20px", "padding-bottom"=>"10px"), children=[
                html_p("local FDR:"),
                # control element for selecting local FDR in current selection for ligation points plots and complementarity plots
                dcc_slider(id="aggregation-fdr-slider", min=0.0, max=1.0, step=0.01, value=1.0,
                    marks=Dict(0.0=>Dict("label"=>"0.0"), 1.0=>Dict("label"=>"1.0")),
                    tooltip=Dict("placement"=>"bottom", "always_visible"=>true)),
            ]
        ),
        html_div(
            className="horizontal", # horizontal is defined in assets/general.css
            children=[
                # Dash component to handle copying data for selected points in ligation points plots and complementarity plots to clipboard
                # this component uses the systems clipboard api, if none is present, it does not work
                dcc_download(id="clip"),
                html_p(id="clip-text", "(click in plot to download tooltip data)"),
            ]
        )

    ]
)

# combine all control components and the info area for additional plots into a column on the left side of the page
control_column_layout(datasets::Vector{String}, min_reads::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64) = html_div(
    id="left-column", # left-column is defined in assets/general.css
    children=[
        headline_layout(),
        html_div(
            className="container", # container is defined in assets/general.css
            children=[
                dataset_and_postioning_layout(datasets),
                reads_selection_layout(min_reads, max_fisher_fdr, max_bp_fdr),
                search_layout(),
                info_area_layout()
            ]
        )
    ]
)

# layout for the cytoscape component
cytoscape_layout(stylesheet::Vector{Dict{String,Any}}) = html_div(
    id="graph-container",
    className="container", # container is defined in assets/general.css
    children=[
        # JavaScript cytoscape component from DashBio loaded as an external Dash component.
        # Documentation available at https://dash.plotly.com/julia/cytoscape/layout
        cytoscape(
            id="graph",
            elements=[],
            autoRefreshLayout=true,
            stylesheet=stylesheet, # load the stylesheet defined in cytostyle.jl
            responsive=true,
            layout=Dict("name"=>"preset", "animate"=>false),
            minZoom=0.2,
            maxZoom=0.9,
            zoom=0.9
        ),
        html_div(
            children=[
                # button for downloading the current network as a .svg file
                html_button("DOWNLOAD SVG", id="save-svg", n_clicks=0),
                html_div(style=Dict("height"=>"8%")),
                # html div for displaying additional info for the network, currently not in use
                html_div(id="legend-container"),
                html_div(style=Dict("height"=>"100%"))
            ]
        )
    ]
)

# internal layout for the circos plot component. Documentation available at https://dash.plotly.com/dash-bio/circos (currently only for python,
# but component properties are the same for julia)
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
# layout for the plotting tab view containing additional plots including the circos plot component
circos_layout(genome_info::Vector{Pair{String,Int}}) = html_div(
    id="plots-container",
    className="container", # container is defined in assets/general.css
    children=[

        html_div(className="horizontal", children=[
            # dropdown to choose the plot type
            dcc_dropdown(id="fdr-source", options=[
                (label = "Node degree distribution", value = "degree"),
                (label = "Annotation stats", value = "annotation"),
                (label = "Odds ratio distribution", value = "odds"),
                (label = "Basepairing alignments clipping distribution", value = "bp"),
            ], value="bp", clearable=false, multi=false),
            html_p("FDR for plots:"),
            # slider to set the FDR cutoff for the plots
            dcc_slider(id="fdr-value", min=0.0, max=1.0, step=0.01, value=0.1,
                marks=Dict(0.0=>Dict("label"=>"0.0"), 1.0=>Dict("label"=>"1.0")),
                tooltip=Dict("placement"=>"bottom", "always_visible"=>true)), #  )
        ], style=Dict("width"=>"100%", "padding-top"=>"20px")),

        # JavaScript circos plot component from DashBio loaded as an external Dash component.
        # Documentation available at https://dash.plotly.com/dash-bio/circos
        html_div(id="circos-container",
        children=[
            html_div(className="plot", children=[dcc_graph(id="plot1")]), # area for hosting the first plot of the selected plot type
            html_div(className="plot", children=[dcc_graph(id="plot2")]), # area for hosting the second plot of the selected plot type
            html_div(className="plot", children=[dcc_graph(id="plot3")]), # area for hosting the general complementarity score distribution plot
            circos(
                id="my-dashbio-circos",
                enableZoomPan=true,
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
        ]),
    ]
)

# layout for the data table component of Dash. Documentation available at https://dash.plotly.com/julia/datatable
table_layout() = html_div(
    id="table-container",
    className="container",
    children=[
        dash_datatable(
            id="table",
            columns=[n == "fdr" ? Dict("name"=>n, "id"=>n, "format"=>Dict("specifier"=>".3f")) : Dict("name"=>n, "id"=>n)
                for n in ["name1", "type1", "name2", "type2", "nb_ints", "fisher_fdr", "odds_ratio", "bp_fdr", "in_libs"]],
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
        # control element for downloading the currently selected interactions as a .csv file
        html_div(children=[
            html_button("DOWNLOAD CSV", id="btn-csv"),
            dcc_download(id="download-dataframe-csv"),
        ], style=Dict("padding-left"=>"20px"))
    ]
)

# layout for the summary view.
summary_layout() = html_div(
    id="summary-container",
    className="container",
    children=[
        html_h3("loading...")  # Defaults to showing a loading message until it is refreshed with actual data
    ]
)

# layout for the how-to view.
howto_layout() = html_div(
    id="howto-container",
    className="container",
    children=[
        html_h3("How to use this interface"),
        html_p("This interactive web application allows you to explore high-throughput RNA-RNA interaction experiments visually. 
            The interface is split into control elements located in the boxes on the left and three main analysis modes located 
            in the tabs on the right. The analysis modes are the graph-based network view, the table view and the plots view. 
            Each of them provides different options to interact with the experimental data."),
        html_table([
            html_tr([
                html_td(className="howto-1", "Control elements"), 
                html_td(className="howto-2", "dataset and node positioning"), 
                html_td(className="howto-3", "Here you can choose the dataset you want to analyze and how the nodes in the network view will be positioned.
                    The clustered positioning employs a stress majorization, which can get slow for selections including more than 500 genes.")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "filtering of the data"), 
                html_td(className="howto-3", "To select a subset of the data, you can set cut-offs for the minimum number of reads an interaction should be found in,
                    the maximum number of interactions to be displayed (interactions are ordered by their read count), and significance levels for the Fisher test of
                    association between two genes (max. fisher fdr) and the comparison of complementarity towards randomly sampled sequnce pairs from the genome of the 
                    studied organism (max. bp fdr). Interactions with no ligation data (paired end reads, where each read maps to a different location but the ligation 
                    site is not included in any of reads) will be displayed with gray edges in the netowrk view, if selected.")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "searching for genes and annotation types"), 
                html_td(className="howto-3", "The selected data can be reduced to interactions involving a single gene or a set of genes. It is also possible to select
                annotation types (as defined in the GFF annotation used in the analysis) instead of genes. The option *exclusive search* defines, how multiple genes or
                annotation types are interpreted. With *exclusive search* enabled, both partners of an interaction have to fit to the selection criteria, while without
                it, only one partner has to fulfil the selected criteria.")
            ]),
            html_tr([
                html_td(className="howto-1", "Network view"), 
                html_td(className="howto-2", "The network"), 
                html_td(className="howto-3", "Here you can adjust the visual representation by moving nodes around and you can select single nodes (genes) or edges
                (interactions) for further information on and plots of the ligation events which were found for the selected entity. The additional plots will be
                displayed on the left below the control elements. For each interaction, a global and a local FDR is computed, the former based on the full dataset,
                the latter on the selected and displayed data.")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "All ligation events between two transcripts"), 
                html_td(className="howto-3", "This plot scatters each ligation point (pair of genomic coordinates of a detected ligation event between two transcripts).
                The plot is interactive and by hovering over each point, additional information is displayed. By clicking into the plot, the underlying data can be
                downloaded as a text file.")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "Ligation summary for a single gene"), 
                html_td(className="howto-3", "This plot aggregates all ligation events involving the selected gene and plots them along the relative coordinates of the
                corresponding genomic region. The plot is interactive and by hovering over it, you can identify the genes ligated to the selected gene at the corresponding
                position. Clicking into it will download all genes found in chimeras at this specific position.")
            ]),
            html_tr([
                html_td(className="howto-1", "Table view"), 
                html_td(className="howto-2", "downloading selected data"), 
                html_td(className="howto-3", "This view summarizes the selected data shown in the graph view as a table. The displayed data can be downloaded using the
                button on the right.")
            ]),
            html_tr([
                html_td(className="howto-1", "Plots view"), 
                html_td(className="howto-2", "Selecting a summary plot"), 
                html_td(className="howto-3", "col3_row2")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "The distribution of complementarity scores"), 
                html_td(className="howto-3", "col3_row1")
            ]),
            html_tr([
                html_td(className="howto-1", ""),
                html_td(className="howto-2", "The circos plot"), 
                html_td(className="howto-3", "col3_row1")
            ]),
            html_tr([
                html_td(className="howto-1", "Summary view"), 
                html_td(className="howto-2", "Total data and selected data"), 
                html_td(className="howto-3", "col3_row2")
            ]),
        ])
    ]
)

# make a tab component from a html div component for integration into a tabs element
astab(div_component::Component, tab_id::String) = dcc_tab(
    id=lowercasefirst(tab_id) * "-tab",
    className="custom-tab",
    selected_className="custom-tab--selected",
    label=uppercasefirst(tab_id),
    value=lowercasefirst(tab_id),
    children=[div_component]
)

# layout of the tabs element hosting the cytoscape graph view, the table view, the additional plots view and the summary view
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
                astab(howto_layout(), "howto"),
            ]
        )
    ]
)

# layout for the complete graphical interface, combining the column with the control elements with the tabs element hosting the
# different data views.
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