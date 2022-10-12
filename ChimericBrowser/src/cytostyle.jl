const core_style = Dict(
    "selector"=>"core",
    "style"=> Dict(
        "selection-box-color"=> "#AAD8FF",
        "selection-box-border-color"=> "#8BB0D0",
        "selection-box-opacity"=> "0.5"
    )
)

const base_node_style = Dict(
    "selector"=> "node",
    "style"=> Dict(
        "width"=> "mapData(current_ratio, 0, 0.7, 55, 150)",
        "height"=> "mapData(current_ratio, 0, 0.7, 55, 150)",
        "content"=> "data(label)",
        "font-size"=> "mapData(current_ratio, 0, 0.7, 12px, 32px)",
        "text-valign"=> "center",
        "text-halign"=> "center",
        "background-color"=> "#555",
        "color"=> "black",
        "overlay-padding"=> "4px",
        "z-index"=> "2",
        "text-wrap"=> "wrap",
        "opacity"=> "0.75"
    )
)

const selected_node_style = Dict(
    "selector"=> "node:selected",
    "style"=> Dict(
        "opacity"=>"1.0",
        "z-index"=>"3",
        "border-width"=> "6px",
        "border-color"=> "#AAD8FF",
        "border-opacity"=> "0.5"
    )
)

const srna_node_style = Dict(
    "selector"=> ".sRNA",
    "style"=> Dict(
        "background-color"=> "CadetBlue",
        "font-weight"=>"bold"
    )
)

const cds_node_style = Dict(
    "selector"=> ".CDS",
    "style"=> Dict(
        "background-color"=> "BurlyWood",
        "font-weight"=>"bold"
    )
)

const utr5_node_style = Dict(
    "selector"=> ".5UTR",
    "style"=> Dict(
        "background-color"=> "BurlyWood",
        "font-weight"=>"bold"
    )
)

const utr3_node_style = Dict(
    "selector"=> ".3UTR",
    "style"=> Dict(
        "background-color"=> "BurlyWood",
        "font-weight"=>"bold"
    )
)

const base_edge_style = Dict(
    "selector"=> "edge",
    "style"=> Dict(
        "curve-style"=> "bezier",
        "opacity"=> "0.4",
        "target-arrow-shape"=> "triangle",
        "line-color"=> "#aaa",
        "width"=> "mapData(current_ratio, 0, 0.7, 2, 20)",
        "z-index"=>"1",
        "overlay-padding"=> "3px"
    )
)

const selected_edge_style = Dict(
    "selector"=> "edge:selected",
    "style"=> Dict(
        "opacity"=>"1.0",
        "z-index"=>"3"
    )
)

const srna_edge_style = Dict(
    "selector"=> ".srna_edge",
    "style"=> Dict(
      "line-color"=> "mapData(relpos, 0, 1, red, blue)"
    )
)

const stylesheet = [core_style, base_node_style, selected_node_style, srna_node_style, cds_node_style, utr5_node_style, utr3_node_style,
    base_edge_style, selected_edge_style, srna_edge_style]