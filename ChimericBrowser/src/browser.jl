function chimeric_browser(results_folder::String, genome_file::String, types::Dict{String,String}, min_reads::Int, max_fdr::Float64, address::String, port::Int)

    interactions, gene_name_info, gene_name_position, genome_info = load_data(results_folder, genome_file, min_reads, max_fdr)

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions)]), genome_info, stylesheet(types))

    update_selection_callback!(app, interactions, gene_name_info, gene_name_position, types["srna"])
    update_dataset_callback!(app, gene_name_info)
    update_selected_element_callback!(app)
    click_cyto_button_callback!(app)
    click_table_button_callback!(app, interactions)
    click_add_node_callback!(app)

    run_server(app, address, port, debug=false)
end
