function chimeric_browser(results_folder::String, genome_file::String, srna_type::String)

    interactions, gene_names_types, gene_name_position, genome_info = load_data(results_folder, genome_file)

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions)]), genome_info)

    update_selection_callback!(app, interactions, gene_names_types, gene_name_position, srna_type)
    update_dataset_callback!(app, gene_names_types)
    update_selected_element_callback!(app)
    click_cyto_button_callback!(app)
    click_table_button_callback!(app, interactions)
    click_add_node_callback!(app)
    #update_layout_callback!(app::Dash.DashApp)

    run_server(app, "0.0.0.0", debug=true)
end