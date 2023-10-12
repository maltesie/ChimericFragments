# Merged alignments with regard to their presence on reads 1 and 2 in paired end illumina sequencing data.
struct MergedAlignedRead
    alns::AlignedReads
    pindexpairs::Vector{Tuple{Int,Int}}
    merged::Vector{Bool}
end

# Empty constructor
MergedAlignedRead(alnreads::AlignedReads) = MergedAlignedRead(alnreads, Tuple{Int, Int}[], Bool[])

# Number of merged alignments in the struct
Base.length(mergedread::MergedAlignedRead) = length(mergedread.pindexpairs)

# Two alignments on one read can be merged, if they are on different reads, have the same annotation and if the leftmost coordinate of
# the alignment on read 1 is upstream or equal to the leftmost coordinate of the alignment on read 2.
canmerge(alns::AlignedReads, i1::Int, i2::Int) = alns.annotated[i2] && (alns.reads[i1]!==alns.reads[i2]) && (alns.annames[i1]===alns.annames[i2]) &&
((alns.strands[i1] == STRAND_POS) ? (alns.leftpos[i1]<=alns.leftpos[i2]) : (alns.rightpos[i1]>=alns.rightpos[i2]))

# Do the merge on preallocated mergedread for efficiency.
function mergeparts!(mergedread::MergedAlignedRead, alnread::AlignedRead)
    alns = alnread.alns

    # Setup a reference to the alignments to be merged
    index = @view alns.pindex[alnread.range]

    # Set preallocated mergedread to appropriate size and overwrite with defaults
    length(mergedread.pindexpairs) >= length(index) || resize!(mergedread.pindexpairs, length(index))
    length(mergedread.merged) >= length(index) || resize!(mergedread.merged, length(index))
    mergedread.merged .= false
    c::Int = 0

    # loop through all alignments
    for (i::Int, ii::Int) in enumerate(index)

        # skip if current alignment was already merged or has no annotation
        (mergedread.merged[i] || !alns.annotated[ii]) && continue

        c += 1
        # Set index of remaining alignments to check for merging
        nextolpi = i+1

        # loop through remaining alignments on read until one can be merged
        if nextolpi <= length(index)
            while !canmerge(alns, ii, index[nextolpi])
                nextolpi += 1
                (nextolpi > length(index)) && break
            end
        end

        # merge if merable alignment was found and mark as merged.
        if nextolpi > length(index)
            mergedread.pindexpairs[c] = (ii, ii)
        else
            mergedread.pindexpairs[c] = (ii, index[nextolpi])
            mergedread.merged[nextolpi] = true
        end
    end
    resize!(mergedread.pindexpairs, c)
    return mergedread
end

# Check the space between alignments on one read to see if there is a ligation point.
function hasligationpoint(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; max_distance=5)
    alns = mergedread.alns
    return ((alns.reads[first(pair1)]==alns.reads[first(pair2)]) && (alns.read_leftpos[first(pair2)]-alns.read_rightpos[first(pair1)] <= max_distance)) ||
    ((alns.reads[last(pair1)]==alns.reads[last(pair2)]) && (alns.read_leftpos[last(pair2)]-alns.read_rightpos[last(pair1)] <= max_distance))
end

# Compute overlap or distance between two intervals, if order is checked, interval 1 has to be upstream of interval 2 for finite distance
function distance(l1::Int, r1::Int, s1::Strand, l2::Int, r2::Int, s2::Strand; check_order=false)::Float64
    s1 != s2 && return Inf
    check_order && (s1 == STRAND_POS ? l1 > l2 : r2 < r1) && return Inf
    l2>r1 && return l2-r1+1
    l1>r2 && return l1-r2+1
    return 0
end

# A pair of alignments is chimeric.
function ischimeric(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; min_distance=1000, check_annotation=true, check_order=false)
    # if enabled, check annotation names
    check_annotation && (mergedread.alns.annames[first(pair1)] == mergedread.alns.annames[first(pair2)]) && return false
    # if enabled, include order of appearance on the reads. If something upstream of something else in the reference appears in the wrong order,
    # the alignments are classified as chimeric. This is how self-chimeras are defined.
    return distance(
        mergedread.alns.leftpos[first(pair1)],
        mergedread.alns.rightpos[last(pair1)],
        mergedread.alns.strands[first(pair1)],
        mergedread.alns.leftpos[first(pair2)],
        mergedread.alns.rightpos[last(pair2)],
        mergedread.alns.strands[first(pair2)]; check_order=check_order) > min_distance
end

# Check if any pair of alignments on a read is chimeric
function ischimeric(mergedread::MergedAlignedRead; min_distance=1000, check_annotation=true, check_order=false)
    length(mergedread) > 1 || return false
    for (i, p1) in enumerate(mergedread.pindexpairs), p2 in (@view mergedread.pindexpairs[i+1:end])
        ischimeric(mergedread, p1, p2; min_distance=min_distance, check_annotation=check_annotation, check_order=check_order) && return true
    end
    return false
end

# Check if there are self-chimeric interactions on the read (see comments in ischimeric)
function isselfchimeric(mergedread::MergedAlignedRead; min_distance=1000)
    length(mergedread) > 1 || return false
    for (i, p1) in enumerate(mergedread.pindexpairs), p2 in (@view mergedread.pindexpairs[i+1:end])
        if !ischimeric(mergedread, p1, p2; min_distance=min_distance, check_annotation=true, check_order=false)
            if ischimeric(mergedread, p1, p2; min_distance=min_distance, check_annotation=false, check_order=true)
                mergedread.alns.annames[p1[1]] == mergedread.alns.annames[p2[1]] && return true
            end
        end
    end
    return false
end

# Check if a read has only one (merged) alignment. If for paired end sequencing only one of two reads map, do not classify as single.
function issingle(mergedread::MergedAlignedRead; is_paired_end=true)
    if (length(mergedread) == 1)
        return is_paired_end != (mergedread.pindexpairs[1][1] == mergedread.pindexpairs[1][2])
    end
    return false
end

# Compute hash from name and type of annotation for efficient indexing
function myhash(alns::AlignedReads, i::Int; use_type=true)
    return use_type ? hash(alns.annames[i], hash(alns.antypes[i])) : hash(alns.annames[i])
end

# Compute length of a merged alignment by combining info from both parts on reads 1 and 2
reflen(alns::AlignedReads, pindexpair::Tuple{Int,Int}) =
    max(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)]) -
    min(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)])