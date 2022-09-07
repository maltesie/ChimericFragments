function demultiplex(testseq::T, barcode_queries::Vector{ApproximateSearchQuery{typeof(isequal), T}}; k=1) where T<:BioSequence
    for i in 1:length(barcode_queries)
        isnothing(findfirst(barcode_queries[i], k, testseq)) || (return i)
    end
    return -1
end
function split_libs(infile1::String, prefixfile::Union{String,Nothing}, infile2::Union{String,Nothing}, libname_to_barcode::Dict{String,K}, output_path::String;
                        bc_len=length(first(values(libname_to_barcode))), check_range=1:bc_len, allowed_barcode_distance=1, overwrite_existing=false) where K<:BioSequence

    barcode_queries = [ApproximateSearchQuery(v) for v in values(libname_to_barcode)]
    for bc in values(libname_to_barcode)
        sum(!isnothing(findfirst(barcode_queries[i], 2*allowed_barcode_distance, bc)) for i in 1:length(libname_to_barcode)) > length(libname_to_barcode) &&
            throw(AssertionError("The supplied barcodes do not support the allowed_barcode_distance!"))
    end

    output_files = isnothing(infile2) ?
    [joinpath(output_path, "$(name).fastq.gz") for name in keys(libname_to_barcode)] :
    [(joinpath(output_path, "$(name)_1.fastq.gz"), joinpath(output_path, "$(name)_2.fastq.gz")) for name in keys(libname_to_barcode)]
    all_present = true
    for file in output_files
        if isnothing(infile2)
            if isfile(file)
                if overwrite_existing
                    rm(file)
                    all_present =false
                end
            else
            	all_present =false
            end
        else
            if isfile(file[1])
                if overwrite_existing
                    rm(file[1])
                    all_present =false
                end
            else
            	all_present =false
            end
            if isfile(file[2])
                if overwrite_existing
                    rm(file[2])
                    all_present =false
                end
            else
            	all_present =false
            end
        end
    end
    all_present && return isnothing(infile2) ? FastqgzFiles(output_files) : PairedSingleTypeFiles(output_files, ".fastq.gz", "_1", "_2")
    nb_stats = length(libname_to_barcode)+1
    stats::Vector{Int} = zeros(Int, nb_stats)
    record1::FASTQ.Record = FASTQ.Record()
    isnothing(infile2) || (record2::FASTQ.Record = FASTQ.Record())
    isnothing(prefixfile) || (recordp::FASTQ.Record = FASTQ.Record())
    endswith(infile1, ".gz") ? reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r"))) : reader1 = FASTQ.Reader(open(infile1, "r"))
    isnothing(infile2) || (endswith(infile2, ".gz") ?
                            reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r"))) :
                            reader2 = FASTQ.Reader(open(infile2, "r")))
    isnothing(prefixfile) || (endswith(prefixfile, ".gz") ?
                                readerp = FASTQ.Reader(GzipDecompressorStream(open(prefixfile, "r"))) :
                                readerp = FASTQ.Reader(open(prefixfile, "r")))
    writers = isnothing(infile2) ?
                [FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2))
                for outfile1 in output_files] :
                [[FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2)),
                FASTQ.Writer(GzipCompressorStream(open(outfile2, "w"), level=2))]
                for (outfile1, outfile2) in output_files]
    isnothing(infile2) ? push!(writers, FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified.fastq.gz"), "w"), level=2))) :
        push!(writers, [FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified_1.fastq.gz"), "w"), level=2)),
                        FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified_2.fastq.gz"), "w"), level=2))])
    c = 0
    while !eof(reader1)
        read!(reader1, record1)
        isnothing(infile2) || read!(reader2, record2)
        isnothing(prefixfile) || read!(readerp, recordp)
        c += 1
        checkseq = isnothing(prefixfile) ? LongDNA{4}(record1.data[record1.sequence])[check_range] : LongDNA{4}(recordp.data[recordp.sequence])
        library_id = demultiplex(checkseq, barcode_queries; k=allowed_barcode_distance)
        library_id == -1 && (library_id = nb_stats)
        stats[library_id] += 1
        isnothing(infile2) ?
        write(writers[library_id], record1) :
        (write(writers[library_id][1], record1); write(writers[library_id][2], record2))
    end
    close(reader1)
    isnothing(infile2) || close(reader2)
    for w in writers
        isnothing(infile2) ?
        close(w) :
        (close(w[1]);close(w[2]))
    end
    libnames = String.(keys(libname_to_barcode))
    sorted_index = sortperm(libnames)
    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libnames[sorted_index], stats[sorted_index])])
    count_string *= "\nnot identifyable - $(stats[end])\n"
    count_string = "Counted $c entries in total:\n\n$count_string\n"
    write(joinpath(dirname(infile1), "demultiplex.log"), count_string)
    return isnothing(infile2) ? FastqgzFiles(output_files) : PairedSingleTypeFiles(output_files, ".fastq.gz", "_1", "_2")
end

function split_libs(infile1::String, infile2::String, libname_to_barcode::Dict{String,K}, output_path::String; bc_len=8, check_range=1:bc_len, overwrite_existing=false) where K<:BioSequence
    split_libs(infile1, nothing, infile2, libname_to_barcode, output_path; bc_len=bc_len, check_range=check_range, overwrite_existing=overwrite_existing)
end

function split_libs(infile::String, libname_to_barcode::Dict{String,K}, output_path::String; bc_len=8, check_range=1:bc_len, overwrite_existing=false) where K<:BioSequence
    split_libs(infile, nothing, nothing, libname_to_barcode, output_path; bc_len=bc_len, check_range=check_range, overwrite_existing=overwrite_existing)
end

"""
    Core dispatch of align_mem, which is a wrapper for bwa-mem2. Takes ´in_file1´ and ´in_file2´ and ´genome_file´ and builds a shell
    command that runs bwa-mem2 with the specified additional parameters and writes the resulting alignments in bam format into ´out_file´
    using samtools.

    # Arguments
    - ´min_score::Int´: defines the minimum score for bwa-mem2 to output alignments for
    - ´match::Int´, ´mismatch::Int´, ´gap_open::Int´, ´gamp_extend::Int´: define the affine gap model that is used for scoring
    - ´clipping_penalty::Int´: scoring penalty for clipping at the beginning or the end of the alignment.
    - ´unpair_penalty::Int´: scoring penalty only used in paired end mapping. Penalty for mappings that do not set the paired flag
    - ´unpair_rescue::Bool´: perform Smith-Waterman to try to rescue read pair if pair is lost.
    - ´min_seed_len::Int´: Minimum seed length. Matches shorter than INT will be missed.
    - ´reseeding_factor::Float64´: Trigger re-seeding for a MEM longer than minSeedLenFLOAT. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy.
"""
function align_mem(in_file1::String, in_file2::Union{String,Nothing}, out_file::String, genome_file::String;
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false,
    min_seed_len=18, reseeding_factor=1.4, is_ont=false, threads=6, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    params = ["-A", match, "-B", mismatch, "-O", gap_open, "-E", gap_extend, "-T", min_score, "-L", clipping_penalty, "-r", reseeding_factor, "-k", min_seed_len, "-t", threads]
    isnothing(in_file2) || append!(params, ["-U", unpair_penalty])
    is_ont && append!(params, ["-x", "ont2d"])
    unpair_rescue && push!(params, "-P")
    fileparams = isnothing(in_file2) ? [genome_file, in_file1] : [genome_file, in_file1, in_file2]
    stats_file = out_file * ".log"
    run(pipeline(
        `$bwa_bin mem -v 1 $params $fileparams`,
        stdout = pipeline(
            `$sam_bin view -u`,
            stdout = pipeline(
                `$sam_bin sort -o $out_file`))))
    run(pipeline(
        `$sam_bin index $out_file`,
        `$sam_bin stats $out_file`,
        stats_file))
end
align_mem(in_file::String, out_file::String, genome_file::String;
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, min_seed_len=18,
    reseeding_factor=1.4, is_ont=false, threads=6, bwa_bin="bwa-mem2", sam_bin="samtools") =
    align_mem(in_file, nothing, out_file::String, genome_file::String;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, clipping_penalty=clipping_penalty, min_seed_len=min_seed_len,
        reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads, bwa_bin=bwa_bin, sam_bin=sam_bin)

"""
    Helper dispatch of align_mem, which is a wrapper for bwa-mem2. Runs align_mem on `read_files::T`
    where T is a FileCollection and aligns it against `genome::Genome`. This enables easy handling of
    whole folders containing sequence files and the manipulation of a genome using Genome.
"""
function align_mem(read_files::T, genome::Genome; min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1,
    clipping_penalty=5, unpair_penalty=9, reseeding_factor=1.4, min_seed_len=18, unpair_rescue=false, is_ont=false,
    threads=6, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:FileCollection}

    tmp_genome = tempname()
    write(tmp_genome, genome)
    outfiles = Vector{String}()
    for file in read_files
        out_file = isa(read_files, SingleTypeFiles) ? file[1:end-length(read_files.type)] * ".bam" : first(file)[1:end-length(read_files.type)-length(read_files.suffix1)] * ".bam"
        push!(outfiles, out_file)
        (isfile(out_file) && !overwrite_existing) && continue
        isa(read_files, SingleTypeFiles) ?
        align_mem(file, out_file, tmp_genome;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open,
                gap_extend=gap_extend, clipping_penalty=clipping_penalty,  min_seed_len=min_seed_len,
                reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(first(file), last(file), out_file, tmp_genome;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
                clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
                min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
                bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        isfile(tmp_genome * ending) && rm(tmp_genome * ending)
    end
    return SingleTypeFiles(outfiles)
end

"""
    Helper dispatch of align_mem, which is a wrapper for bwa-mem2. Runs align_mem on `reads::T` where T
    is a SequenceContainer and aligns it against `genomes::Vector{Genome}`. This enables easy handling
    of Sequences and their alignment against multiple genomes.
"""
function align_mem(reads::Sequences, genomes::Vector{Genome}, out_file::String; min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1,
    clipping_penalty=5, unpair_penalty=9, unpair_rescue=false, min_seed_len=18, reseeding_factor=1.4, is_ont=false, threads=6,
    bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false)

    (isfile(out_file) && !overwrite_existing) && return
    tmp_reads = tempname()
    tmp_reads2 = tempname()
    tmp_genome = tempname()
    ispaired = reads.seqnames[1:2:end] == reads.seqnames[2:2:end]
    if ispaired
        write(tmp_reads, tmp_reads2, reads)
    else
        write(tmp_reads, reads)
    end
    for (i,genome) in enumerate(genomes)
        write(tmp_genome, genome)
        this_out_file = out_file
        length(genomes) > 1 && (this_out_file = joinpath(dirname(out_file), "$(i)_" * basename(out_file)))

        !ispaired ?
        align_mem(tmp_reads, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor,
            is_ont=is_ont, threads=threads, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(tmp_reads, tmp_reads2, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
            min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
            bwa_bin=bwa_bin, sam_bin=sam_bin)

        rm(tmp_genome)
    end
    rm(tmp_reads)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
    if length(genomes) > 1
        bam_files = []
        bai_files = []
        for i in 1:length(genomes)
            push!(bam_files, joinpath(dirname(out_file), "$(i)_" * basename(out_file)))
            push!(bai_files, joinpath(dirname(out_file), "$(i)_" * basename(out_file) * ".bai"))
        end
        run(`$sam_bin merge -X $out_file $bam_files $bai_files`)
        run(`$sam_bin index $out_file`)
        for (bam_file, bai_file) in zip(bam_files, bai_files)
            rm(bam_file)
            rm(bai_file)
        end
    end
end
align_mem(reads::T, genome::Genome, out_file::String;
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, unpair_penalty=9, min_seed_len=18, reseeding_factor=1.4,
    unpair_rescue=false, is_ont=false, threads=6, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:Sequences} =
    align_mem(reads, [genome], out_file;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, unpair_penalty=unpair_penalty,
        unpair_rescue=unpair_rescue, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
        bwa_bin=bwa_bin, sam_bin=sam_bin, overwrite_existing=overwrite_existing)
