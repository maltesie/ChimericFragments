function Genome(sequence_dict::Dict{String, LongSequence})
    seq = LongDNA{4}(undef, sum(length(s) for s in values(sequence_dict)))
    chrs = Dict{String, UnitRange}()
    si = 1
    for (name, sequence) in sequence_dict
        srange = si:si+length(sequence)-1
        seq[srange] = sequence
        push!(chrs, name=>srange)
        si = last(srange) + 1
    end
    return Genome(seq, chrs)
end
Genome(sequences::Vector{LongSequence}, names::Vector{String}) = Genome(Dict(n=>s for (n,s) in zip(sequences, names)))
Genome(sequence::LongSequence, name::String) = Genome([sequence], [name])
Genome(genome_file::String) = Genome(read_genomic_fasta(genome_file))

Base.length(genome::Genome) = length(genome.seq)
Base.getindex(genome::Genome, key::String) = genome.seq[genome.chroms[key]]
Base.getindex(genome::Genome, key::Pair{String, UnitRange}) = genome[first(key)][last(key)]

function chomosomecount(genome::Genome)
    return length(genome.chroms)
end

function Base.iterate(genome::Genome)
    (chr, slice) = first(genome.chroms)
    ((chr, genome.seq[slice]), 1)
end

function Base.iterate(genome::Genome, state::Int)
    state += 1
    state > genome.chroms.count && (return nothing)
    for (i, (chr, slice)) in enumerate(genome.chroms)
        (i == state) && (return ((chr, genome.seq[slice]), state))
    end
end

function Base.merge(genomes::Vector{Genome})
    length(genomes) == 1 && return genomes[1]
    merged = genomes[1]
    for genome in genomes[2:end]
        merged *= genome
    end
    return merged
end

function Base.:*(genome1::Genome, genome2::Genome)
    return Genome(genome1.seq*genome2.seq, merge(genome1.chroms, Dict(key=>(range .+ length(genome1)) for (key, range) in genome2.chroms)))
end

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>seq for (chr, seq) in genome), file)
end

function read_genomic_fasta(fasta_file::String)
    genome::Dict{String, LongSequence} = Dict()
    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            genome[identifier(record)] = FASTA.sequence(record)
        end
    end
    return genome
end

function write_genomic_fasta(genome::Dict{String, T}, fasta_file::String) where T <: BioSequence
    open(FASTA.Writer, fasta_file) do writer
        for (chr, seq) in genome
            write(writer, FASTA.Record(chr, seq))
        end
    end
end
