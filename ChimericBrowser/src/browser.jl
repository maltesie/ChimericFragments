function chimeric_browser(results_folder::String, genome_file::String, srna_type::String)

    interactions_dfs, singles_dfs, stats_dfs, gene_names_types, genome_info = load_data(results_folder, genome_file)

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions_dfs)]), [Dict("label"=>"test", "value"=>"test")], genome_info)

    update_selection_callback!(app, interactions_dfs, srna_type, gene_names_types)
    update_dataset_callback!(app, gene_names_types)
    update_selected_element_callback!(app, stats_dfs)

    run_server(app, "0.0.0.0", debug=true)
end