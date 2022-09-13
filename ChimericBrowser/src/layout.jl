cyto_layout = html_div([
    cytoscape(
        id="cytoscape-two-nodes",
        layout=Dict("name" =>  "preset"),
        style=Dict("width" =>  "100%", "height" =>  "400px"),
        elements=[
            Dict("data" =>  Dict("id" =>  "one", "label" =>  "Node 1"), "position" =>  ("x" =>  75, "y" =>  75)),
            Dict("data" =>  Dict("id" =>  "two", "label" =>  "Node 2"), "position" =>  ("x" =>  200, "y" =>  200)),
            Dict("data" =>  Dict("source" =>  "one", "target" =>  "two"))
        ]
    )
])

circos_layout = html_div([
    circos(
        id="circos-simple",
        layout=[
            Dict(
                "id"=> "NC_002505",
                "label"=> "",
                "color"=> "#009933",
                "len"=> 2961149
            ),
            Dict(
                "id"=> "NC_002506",
                "label"=> "",
                "color"=> "#99cc33",
                "len"=> 1072315
            )
        ]
    )
])

browser_layout = html_div([circos_layout, cyto_layout])