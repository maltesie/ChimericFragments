function chimeric_browser(results_folder::String, genome_file::String, types::Dict{String,String}, min_reads::Int, max_fdr::Float64, max_bp_fdr::Float64,
        address::String, port::Int, check_interaction_distances::Tuple{Int,Int}, param_dict::Vector{Pair{String,String}}, bp_parameters::NTuple{6,Int})

    interactions, genome_info, genome = load_data(results_folder, genome_file, min_reads, max_fdr, max_bp_fdr)

    scores = Dict(
        (DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
        (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
        (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3]
    )
    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
                                                            gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions)]), genome_info, stylesheet(types), min_reads, max_fdr, max_bp_fdr)

    update_selection_callback!(app, interactions, types["srna"], param_dict)
    update_dataset_callback!(app, interactions, min_reads)
    update_selected_element_callback!(app, genome, interactions, check_interaction_distances, bp_parameters)
    click_cyto_button_callback!(app)
    click_table_button_callback!(app, interactions)
    click_add_node_callback!(app)

    run_server(app, address, port, debug=false)
end
