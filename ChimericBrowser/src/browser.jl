# Compute the complementarity score as a combination between the alignment score from the AffineGapScoreModel and the weighted shift between
# complementarity regions towards eachother in the local coordinates of the fixed-length sequences used to compute the complementarity.
score_bp(paln::PairwiseAlignmentResult, shift_weight::Float64) = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))

# Main function for building the visualization. Loads and initializes data and sets up Dash app with layout and callbacks.
function chimeric_browser(results_folder::String, genome_file::String, types::Dict{String,String}, min_reads::Int, max_fisher_fdr::Float64, max_bp_fdr::Float64,
        address::String, port::Int, check_interaction_distances::Tuple{Int,Int}, param_dict::Vector{Pair{String,String}}, bp_parameters::NTuple{6,Int},
        n_genome_samples::Int, shift_weight::Float64)

    # Load Interactions struct from jld2 files generated by analyze.jl and read genome.
    interactions, genome_info, genome = load_data(results_folder, genome_file, min_reads, max_fisher_fdr, max_bp_fdr)

    # Setup AffineGapScoreModel for computing complementarity between sequences. Parameters are defined in config file passed to visualize.jl
    scores = Dict(
        (DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
        (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
        (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3]
    )
    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
                                                            gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

    # Compute length from interval around ligation points. Sequences are taken from the genome.
    seq_length = check_interaction_distances[1]-check_interaction_distances[2]

    # Precompute complement and reversed genomes
    genomeseq = foldl(*, values(genome))
    cgenomeseq = BioSequences.complement(genomeseq)
    rgenomeseq = BioSequences.reverse(genomeseq)
    rcgenomeseq = BioSequences.complement(rgenomeseq)

    # Compute empirical density of scores from pairs of fixed-length random sequences from the genome to plot random model used for computing p-values
    genome_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
        i1 % 2 == 0 ? view(genomeseq, i1:i1+seq_length-1) : view(rcgenomeseq, i1:i1+seq_length-1),
        i2 % 2 == 0 ? view(rgenomeseq, i2:i2+seq_length-1) : view(cgenomeseq, i2:i2+seq_length-1),
        model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genomeseq)-seq_length), (n_genome_samples,2)))]
    )

    # Compute empirical density of scores from pairs of random sequences to plot expectation for totally random sequence pairs.
    randseq = BioSequences.randdnaseq(length(genomeseq))
    randseq_model_ecdf = ecdf([score_bp(pairalign(LocalAlignment(),
        view(randseq, i1:i1+seq_length-1),
        view(randseq, i2:i2+seq_length-1),
        model), shift_weight) for (i1, i2) in eachrow(rand(1:(length(genomeseq)-seq_length), (n_genome_samples,2)))]
    )

    # Setup Dash app and connect with HTML layout defined in layout.jl
    app = dash(assets_folder=assets_path)
    app.layout = browser_layout(sort([k for k in keys(interactions)]), genome_info, stylesheet(types), min_reads, max_fisher_fdr, max_bp_fdr)

    # Connect callback functions defined in callbacks.jl
    update_selection_callback!(app, interactions, param_dict)
    update_dataset_callback!(app, interactions, min_reads)
    update_plots_callback!(app, interactions, randseq_model_ecdf, genome_model_ecdf)
    update_aggregation_slider_callback!(app, genome, interactions, check_interaction_distances, model)
    update_selected_element_callback!(app)
    click_clipboard_callback!(app, interactions, check_interaction_distances[1])
    click_cyto_button_callback!(app)
    click_table_button_callback!(app, interactions)
    click_add_node_callback!(app)

    # Start webserver hosting the Dash app
    run_server(app, address, port, debug=false)
end
