type(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.type
name(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.name
refname(feature::Interval) = feature.seqname
Base.length(i::Interval{T}) where {T<:AnnotationStyle} = rightposition(i) - leftposition(i) + 1
featureparams(feature::Interval{Annotation}) = feature.metadata.params
featureparam(feature::Interval{Annotation}, key::String) = feature.metadata.params[key]
featureparam(feature::Interval{Annotation}, key::String, ::Type{T}) where {T} = parse(T, feature.metadata.params[key])
setfeatureparam!(feature::Interval{Annotation}, key::String, value::String) = feature.metadata.params[key] = value
hasannotationkey(feature::Interval{Annotation}, key::String) = key in keys(featureparams(feature))
annotation(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata
hasannotation(feature::Interval{T}) where T<:AnnotationStyle = !isempty(annotation(feature))
ispositivestrand(feature::Interval{T}) where T<:AnnotationStyle = strand(feature) === STRAND_POS

function summarize(feature::Interval{Annotation})
    s = "Feature: [$(leftposition(feature)), $(rightposition(feature))] on $(refname(feature)) ($(strand(feature))) - "
    s *= hasannotation(feature) ? "annotation: $(type(feature)):$(name(feature))" : "not annotated"
    s *= "\nparameters: $(paramstring(feature))\n"
end
function Base.show(feature::Interval{Annotation})
    println(summarize(feature))
end

function Annotation()
    Annotation("", "", Dict{String, String}())
end

function Annotation(type::String, name::String; args...)
    Annotation(type, name, merge(Dict("Name"=>name), Dict{String, String}("$key"=>value for (key,value) in args)))
end

Base.isempty(annotation::Annotation) = isempty(annotation.type) && isempty(annotation.name)

function AlignmentAnnotation()
    AlignmentAnnotation("", "", 0)
end

Base.isempty(annotation::AlignmentAnnotation) = isempty(annotation.type) && isempty(annotation.name) && (annotation.overlap==0)

name(annot::T) where {T<:AnnotationStyle} = annot.name
type(annot::T) where {T<:AnnotationStyle} = annot.type
overlap(annot::AlignmentAnnotation) = annot.overlap

function BaseAnnotation(feature::Interval{Annotation}, base_coverage::BaseCoverage)
    left = leftposition(feature)
    right = rightposition(feature)
    seq = base_coverage.ref_seq[left:right]
    ispositivestrand(feature) || reverse_complement!(seq)
    count = ispositivestrand(feature) ? errorcov.fcount : errorcov.rcount
    r = ispositivestrand(feature) ? (left:right) : (right:-1:left)
    ref = Int[(seq[i] in (DNA_A, DNA_T, DNA_G, DNA_C)) ? count[seq[i]][ii] : 0 for (i, ii) in enumerate(r)]
    BaseAnnotation(type(feature), name(feature), ref, count[ref][:A][r], count[:T][r], count[:G][r], count[:C][r], count[:Gap][r], count[:Ins][r])
end

Features() = Features(Annotation)

function Features(::Type{T}) where T
    return Features(Vector{Interval{T}}())
end

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true), Dict{String, Int}())
end

function Features(feature_list::Vector{Interval{Annotation}}, chroms::Dict{String, Int})
    return Features(IntervalCollection(feature_list, true), chroms)
end

function Features(gff_file::String, type::Vector{String}; name_keys=["Name"], same_name_rule=:all)
    same_name_rule in (:first, :all, :none) || throw(AssertionError("same_name_rule must be :first, :all, or :none"))
    features = GFF3.Reader(open(gff_file), skip_directives=false)
    intervals = Vector{Interval{Annotation}}()
    names = Dict{Tuple{String,String},Int}()
    chroms = Dict{String, Int}()
    for feature in features
        if GFF3.isdirective(feature)
            line = GFF3.content(feature)
            if startswith(line, "sequence-region")
                seqid, _, rp = split(line)[2:4]
                push!(chroms, seqid => parse(Int, rp))
            end
            continue
        end
        GFF3.featuretype(feature) == "Source" && continue
        (GFF3.featuretype(feature) in type || isempty(type)) || continue
        seqn = GFF3.seqid(feature)
        name = ("NA", "NA")
        as = GFF3.attributes(feature)
        for name_key in reverse(name_keys)
            for (key, values) in as
                name_key === key || continue
                name = (GFF3.featuretype(feature), join(values, ","))
            end
        end
        name in keys(names) ? (names[name]+=1) : (names[name]=1)
        (same_name_rule === :first && names[name] > 1) && continue
        n = names[name] > 1 ? name[2] * "$(names[name])" : name[2]
        annot = Annotation(GFF3.featuretype(feature), n, Dict(pair[1] => join(pair[2], ",") for pair in GFF3.attributes(feature)))
        push!(intervals, Interval(seqn, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annot))
    end
    same_name_rule === :none && (return Features([i for i in intervals if ((type(i), name(i)) in names && names[(type(i), name(i))] == 1)]))
    return Features(intervals, chroms)
end

function Features(gff_file::String, type::String; name_keys=["Name"], same_name_rule=:all)
    return Features(gff_file, [type], name_keys=name_keys, same_name_rule=same_name_rule)
end

function Features(gff_file::String; name_keys=["Name"], same_name_rule=:all)
    return Features(gff_file, String[], name_keys=name_keys, same_name_rule=same_name_rule)
end

function Features(gff_file::String, bam_file::String, genome::Genome)
    new_intervals = Interval{BaseAnnotation}[]
    features = Features(gff_file)
    errorcov = BaseCoverage(bam_file, genome)
    for feature in features
        annotation = BaseAnnotation(feature, errorcov)
        push!(new_intervals, Interval(refname(feature), leftposition(feature), rightposition(feature), strand(feature), annotation))
    end
    return Features(IntervalCollection(new_intervals, true), genome.chroms)
end

types(features::Features) = Set(type(f) for f in features)
refnames(features::Features) = collect(keys(features.list.trees))
function isoverlapping(alignmentinterval::Interval{T}, feature::Interval{K}) where {T<:AnnotationStyle, K<:AnnotationStyle}
    refname(feature) == refname(alignmentinterval) || return false
    strand(feature) == strand(alignmentinterval) || return false
    return leftposition(feature) <= rightposition(alignmentinterval) && leftposition(alignmentinterval) <= rightposition(feature)
end

Base.push!(features::Features, interval::Interval) = push!(features.list, interval)
function Base.merge!(features1::Features{T}, features2::Features{T}) where T
    for feature in features2
        push!(features1, feature)
    end
end
function Base.merge(features1::Features{T}, features2::Features{T}) where T
    re = Features(T)
    for feature in features1
        push!(re, feature)
    end
    for feature in features2
        push!(re, feature)
    end
    return re
end
Base.:*(features1::Features, features2::Features) = merge(features1, features2)

function mergetypes(features::Features, types::Vector{String}, mergetype::String)
    merged_features = Interval{Annotation}[]
    feature_collector = Dict{String,Vector{Interval{Annotation}}}(name(feature)=>Interval{Annotation}[] for feature in features if type(feature) in types)
    for feature in features
        if type(feature) in types
            push!(feature_collector[name(feature)], feature)
        else
            push!(merged_features, feature)
        end
    end
    for (n, fs) in feature_collector
        any((refname(f)!=refname(fs[1])) || (strand(f)!=strand(fs[1])) for f in fs) &&
            throw(AssertionError("Features with the same name of types $types have to be on the same reference sequence and strand."))
        left = minimum(leftposition(f) for f in fs)
        right = maximum(rightposition(f) for f in fs)
        push!(merged_features, Interval(refname(fs[1]), left, right, strand(fs[1]), Annotation(mergetype, n)))
    end
    return Features(merged_features)
end

Base.iterate(features::T) where T<:AnnotationContainer = iterate(features.list)
Base.iterate(features::T, state::Tuple{Int64,GenomicFeatures.ICTree{I},GenomicFeatures.ICTreeIteratorState{I}}) where {T<:AnnotationContainer,I<:AnnotationStyle} = iterate(features.list, state)
Base.length(features::T) where T<:AnnotationContainer = length(features.list)
Base.split(features::T) where T<:AnnotationContainer = [Features([feature for feature in features if type(feature)==t], [t]) for t in types(features)]

function Base.getindex(features::T, i::Int) where T<:AnnotationContainer
    (i < 1 || i > length(features)) && throw(AssertionError("Trying to access $(length(features))-element $T at index $i"))
    for (c, f) in enumerate(features)
        c === i && return f
    end
end

Base.convert(::Type{Interval{Annotation}}, i::Interval{Nothing}) = Interval(refname(i), leftposition(i), rightposition(i), strand(i), Annotation())
strand_filter(a::Interval, b::Interval)::Bool = strand(a) === strand(b)
function GenomicFeatures.eachoverlap(feature::Interval{T}, features::Features{Annotation}) where T
    T !== Annotation && (feature = Interval{Annotation}(feature))
    haskey(features.list.trees, refname(feature)) ?
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter,
                GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refname(feature)], feature)) :
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter,
                GenomicFeatures.ICTreeIntersection{Annotation}(), GenomicFeatures.ICTree{Annotation}(), feature))
end

function hasoverlap(feature::Interval, features::Features)
    for _ in eachoverlap(feature, features)
        return true
    end
    return false
end
function firstoverlap(feature::Interval, features::Features)
    for int in eachoverlap(feature, features)
        return int
    end
    return nothing
end
function lastoverlap(feature::Interval, features::Features)
    olpinterval = copy(feature)
    for int in eachoverlap(feature, features)
        olpinterval = int
    end
    return olpinterval
end

function paramstring(params::Dict{String,String}; priority=("Name",))
    ps = join(("$key=$(params[key])" for key in priority if key in keys(params)), "; ")
    os = join(("$key=$(params[key])" for key in sort(collect(keys(params))) if !(key in priority)), "; ")
    both = ps * ((isempty(ps) || isempty(os)) ? "" : "; ") * os
    return isempty(both) ? "none" : both
end
paramstring(feature::Interval{Annotation}; priority=("Name",)) = paramstring(featureparams(feature); priority=priority)

function Base.write(file::String, features::Features; zip=false, tabix=false)
    chroms = copy(features.chroms)
    writer = GFF3.Writer(open(file, "w"))
    write(writer, GFF3.Record("##gff-version 3.2.1"))
    for feature in features
        #println(feature)
        if refname(feature) in keys(chroms)
            write(writer, GFF3.Record("##sequence-region $(refname(feature)) 1 $(chroms[refname(feature)])"))
            delete!(chroms, refname(feature))
        end
        write(writer, GFF3.Record("$(refname(feature))\t.\t$(type(feature))\t$(feature.first)\t$(feature.last)\t.\t$(feature.strand)\t.\t$(paramstring(featureparams(feature)))"))
    end
    b = close(writer)
    sleep(0.5)
    zip && run(`bgzip --force $file`)
    tabix && run(`tabix $(file).gz`)
    return b
end

utrright(f::Interval{Annotation}, features::Features, len::Int, utr_type::String, cds_type::String) =
    Interval(refname(f), rightposition(f)+1,
        hasoverlap(Interval(refname(f), rightposition(f)+1, rightposition(f)+2*len, strand(f)), features) ?
            minimum(type(i) == cds_type ? floor(Int, (leftposition(i) + rightposition(f))/2) : min(rightposition(f)+len, leftposition(i)-1)
                for i in eachoverlap(Interval(refname(f), rightposition(f)+1, rightposition(f)+2*len, strand(f)), features)) :
            rightposition(f)+len,
        strand(f), Annotation(utr_type, name(f)))

utrleft(f::Interval{Annotation}, features::Features, len::Int, utr_type::String, cds_type::String) =
    Interval(refname(f),
        hasoverlap(Interval(refname(f), leftposition(f)-2*len, leftposition(f)-1, strand(f)), features) ?
            maximum(type(i) == cds_type ? floor(Int, (leftposition(f) + rightposition(i))/2)+1 : max(leftposition(f)-len, rightposition(i)+1)
                for i in eachoverlap(Interval(refname(f), leftposition(f)-2*len, leftposition(f)-1, strand(f)), features)) :
            leftposition(f)-len,
        leftposition(f)-1, strand(f), Annotation(utr_type, name(f)))

function add5utrs!(features::Features; cds_type="CDS", utr_type="5UTR", utr_length=200, min_utr_length=20)
    new_features = Vector{Interval{Annotation}}()
    for feature in features
        type(feature) != cds_type && continue
        new_feature = ispositivestrand(feature) ?
            utrleft(feature, features, utr_length, utr_type, cds_type) :
            utrright(feature, features, utr_length, utr_type, cds_type)
            ((length(new_feature) < min_utr_length) || (leftposition(new_feature) < 1)) && continue
        push!(new_features, new_feature)
    end
    merge!(features, Features(new_features))
end

function add3utrs!(features::Features; cds_type="CDS", utr_type="3UTR", utr_length=200, min_utr_length=20)
    new_features = Vector{Interval{Annotation}}()
    for feature in features
        type(feature) != cds_type && continue
        new_feature = ispositivestrand(feature) ?
            utrright(feature, features, utr_length, utr_type, cds_type) :
            utrleft(feature, features, utr_length, utr_type, cds_type)
        ((length(new_feature) < min_utr_length) || (leftposition(new_feature) < 1)) && continue
        push!(new_features, new_feature)
    end
    merge!(features, Features(new_features))
end

function addutrs!(features::Features; cds_type="CDS", five_type="5UTR", three_type="3UTR", five_length=200, three_length=200, min_length=20)
    add5utrs!(features; cds_type=cds_type, utr_type=five_type, utr_length=five_length, min_utr_length=min_length)
    add3utrs!(features; cds_type=cds_type, utr_type=three_type, utr_length=three_length, min_utr_length=min_length)
end

function addigrs!(features::Features; igr_type="IGR", min_igr_length=20)
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (feature.strand == STRAND_NEG)]
    for base_features in (base_features_pos, base_features_neg)
        nb_features = base_features === base_features_pos ? length(base_features_pos) : length(base_features_neg)
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            refname(feature) == refname(next_feature) || continue
            stop, start = leftposition(next_feature), rightposition(feature)
            (stop-1) - (start + 1) > min_igr_length || continue
            igr = Interval(refname(feature), start+1, stop-1, base_features === base_features_pos ? STRAND_POS : STRAND_NEG,
                Annotation(igr_type, name(feature)*":"*name(next_feature)))
            push!(new_features, igr)
        end
    end
    merge!(features, Features(new_features))
end

function summarize(features::Features)
    chrs = refnames(features)
    ts = sort(collect(types(features)))
    stats = Dict(chr=>Dict(t=>0 for t in ts) for chr in chrs)
    nb_pos = 0
    nb_neg = 0
    for f in features
        stats[refname(f)][type(f)] += 1
        nb_pos += strand(f) === STRAND_POS
        nb_neg += strand(f) === STRAND_NEG
    end
    printstring = "$(nb_pos+nb_neg) features in total with $nb_pos on + strand and $nb_neg on - strand:\n"
    infotable = DataFrame("type"=>ts, [chr=>[stats[chr][t] for t in ts] for chr in chrs]...)
    printstring *= DataFrames.pretty_table(String, infotable, nosubheader=true)
    return printstring
end

function Base.show(features::Features)
    println(summarize(features))
end
