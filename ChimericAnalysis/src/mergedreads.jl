struct MergedAlignedRead
    alns::AlignedReads{UInt}
    pindexpairs::Vector{Tuple{Int,Int}}
    merged::Vector{Bool}
end

MergedAlignedRead(alnreads::AlignedReads{UInt}) = MergedAlignedRead(alnreads, Tuple{Int, Int}[], Bool[])

Base.length(mergedread::MergedAlignedRead) = length(mergedread.pindexpairs)

canmerge(alns::AlignedReads, i1::Int, i2::Int) = alns.annotated[i2] && (alns.reads[i1]!==alns.reads[i2]) && (alns.annames[i1]===alns.annames[i2]) &&
((alns.strands[i1] == STRAND_POS) ? (alns.leftpos[i1]<=alns.leftpos[i2]) : (alns.rightpos[i1]>=alns.rightpos[i2]))

function mergeparts!(mergedread::MergedAlignedRead, alnread::AlignedRead)
    alns = alnread.alns
    index = @view alns.pindex[alnread.range]
    length(mergedread.pindexpairs) >= length(index) || resize!(mergedread.pindexpairs, length(index))
    length(mergedread.merged) >= length(index) || resize!(mergedread.merged, length(index))
    mergedread.merged .= false
    c::Int = 0
    for (i::Int, ii::Int) in enumerate(index)
        (mergedread.merged[i] || !alns.annotated[ii]) && continue
        c += 1
        nextolpi = i+1

        if nextolpi <= length(index)
            while !canmerge(alns, ii, index[nextolpi])
                nextolpi += 1
                (nextolpi > length(index)) && break
            end
        end

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

function hasligationpoint(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; max_distance=5)
    alns = mergedread.alns
    return ((alns.reads[first(pair1)]==alns.reads[first(pair2)]) && (alns.read_leftpos[first(pair2)]-alns.read_rightpos[first(pair1)] <= max_distance)) ||
    ((alns.reads[last(pair1)]==alns.reads[last(pair2)]) && (alns.read_leftpos[last(pair2)]-alns.read_rightpos[last(pair1)] <= max_distance))
end

function distance(l1::Int, r1::Int, s1::Strand, l2::Int, r2::Int, s2::Strand; check_order=false)::Float64
    s1 != s2 && return Inf
    check_order && (s1 == STRAND_POS ? l1 > l2 : r2 < r1) && return Inf
    l2>r1 && return l2-r1+1
    l1>r2 && return l1-r2+1
    return 0
end

function ischimeric(mergedread::MergedAlignedRead, pair1::Tuple{Int,Int}, pair2::Tuple{Int, Int}; min_distance=1000, check_annotation=true, check_order=false)
    check_annotation && (mergedread.alns.annames[first(pair1)] == mergedread.alns.annames[first(pair2)]) && return false
    return distance(
        mergedread.alns.leftpos[first(pair1)],
        mergedread.alns.rightpos[last(pair1)],
        mergedread.alns.strands[first(pair1)],
        mergedread.alns.leftpos[first(pair2)],
        mergedread.alns.rightpos[last(pair2)],
        mergedread.alns.strands[first(pair2)]; check_order=check_order) > min_distance
end

function ischimeric(mergedread::MergedAlignedRead; min_distance=1000, check_annotation=true, check_order=false)
    length(mergedread) > 1 || return false
    for (i, p1) in enumerate(mergedread.pindexpairs), p2 in (@view mergedread.pindexpairs[i+1:end])
        ischimeric(mergedread, p1, p2; min_distance=min_distance, check_annotation=check_annotation, check_order=check_order) && return true
    end
    return false
end

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

function issingle(mergedread::MergedAlignedRead; is_paired_end=true)
    if (length(mergedread) == 1)
        return is_paired_end != (mergedread.pindexpairs[1][1] == mergedread.pindexpairs[1][2])
    end
    return false
end

function myhash(alns::AlignedReads, i::Int; use_type=true)
    return use_type ? hash(alns.annames[i], hash(alns.antypes[i])) : hash(alns.annames[i])
end

reflen(alns::AlignedReads, pindexpair::Tuple{Int,Int}) =
    max(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)]) -
    min(alns.leftpos[first(pindexpair)], alns.leftpos[last(pindexpair)], alns.rightpos[first(pindexpair)], alns.rightpos[last(pindexpair)])