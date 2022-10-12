function chimeric_browser(results_folder::String, genome_file::String)

    interactions, gene_names_types, genome_info = load_data(results_folder, genome_file)

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions)]), [Dict("label"=>"test", "value"=>"test")], genome_info)

    update_selection_callback!(app, interactions, gene_names_types)
    update_dataset_callback!(app, gene_names_types)
    update_selected_element_callback!(app)

    run_server(app, "0.0.0.0", debug=true)
end