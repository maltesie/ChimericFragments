# FILE PARAMETERS

# If BWA or samtools are not installed systemwide, please provide full paths to the binaries
bwa_mem2_bin = nothing
samtools_bin = nothing
fastp_bin = nothing

# If the annotation file, the genome file and the data folder containing the sample files are not
# located next to the config file, please provide full paths.
annotation_file = "annotation.gff"
genome_file = "genome.fa"
data_folder = "data"

# BWA can handle single end, paired end (2 files) and interleaved paired end (1 file) sequence
# files. For paired end and interleaved paired end files, the names of the first and the second
# read have to be identical.
#
# WARNING: Interleaved sequence files from NCBI have different names for first and second read as
# of November, 2022 and are not supported by bwa-mem2 right now.
is_paired_end = true
is_interleaved_paired_end = false

# If there are two files per sample, here you can set the suffixes that differentiate the files.
# The suffixes are ignored if is_paired_end is false or if is_interleaved_paired_end is true.
suffix_read1 = "_1"
suffix_read2 = "_2"

# All sequence files used for one analysis run have to have the same file type which you can set here.
file_type = ".fastq.gz"

# Here, file names have to be related to conditions. Together with the suffixes and the file type,
# samplename_condition defines, which files are used in this analysis run. As an example, with the
# default settings, ChimericFragments would look for the files
#
# "conditionA_rep1_1.fastq.gz" and "conditionaA_rep1_2.fastq.gz"
#
# in the folder "data" next to this config file and assign them to the condition "conditionA" and
# do the same for all sample names listed under samplename_condition. If is_paired_end is false or
# if is_interleaved_paired_end is true, ChimericFragments will look for one file per sample name,
# ingnoring the suffices (leading to "conditionA_rep1.fastq.gz" beeing assigned to "conditionA"
# in the default settings)
samplename_condition = [
    ("hfq_lcd_1", "hfq_lcd"),
    ("hfq_lcd_2", "hfq_lcd"),
    ("hfq_hcd_1", "hfq_hcd"),
    ("hfq_hcd_2", "hfq_hcd")
]


# FASTP PARAMETERS

# indicate if trimming should be performed
skip_trimming = false

# fastp filters duplicated reads based on the distribution of exact copies.
deduplicate = false

# Cutoff for the average quality of a sliding window that decides if reads are trimmed at the front
# and at the end. Reads shorter than min_length after trimming will be discarded.
average_window_quality=20
min_length=20


# BWA-MEM PARAMETERS

# bwa-mem seeding parameter. No alignments shorter than min_seed_length can be found, can effect
# performance if set to lower values. reseeding_factor can be used to tune accuracy vs performance
# with lower values leading to higher accuracy.
min_seed_len=18
reseeding_factor=1.4

# Smith-Waterman extension paramerter. Only alignments with scores greater or equal to min_score are saved.
# The score of an alignment is the number of matches times match_score minus the number of mismatches, gaps,
# the length of the gaps and clipping at both ends times the respective penalty.
min_alignment_score=20
match_score=1
mismatch_penalty=4
gap_open_penalty=6
gap_extend_penalty=1
clipping_penalty=5

# unpair_penalty and unpair_rescue determine if PE reads get paired or not. If the penalty of the insert
# between the reads is higher than unpair_penalty, the reads get unpaired. unpair_rescue determines, if
# pairing should only be attempted, if one read has no hit or always.
unpair_penalty=9
unpair_rescue=true

# set true, if the alignment files need to be sorted and indexed with a .bai for further processing.
# Does not affect this analysis.
sort_and_index_bam=false

# set number of threads for bwa-mem2
threads=8


# CHIMERIC FRAGMENTS PARAMETERS

# The minimum distance between two alignments on the same strand to be called chimeric. The distance
# computation is aware of the order of the alignments on the insert. I.e., if for alignments A and B,
# A sits further towards the 5' end than B on the insert but their origin locations are reversed, the
# distance between them is set to infinity.
min_distance=1000

# Set true, if the in formation contained in the reads is reverse complementary to its origin on
# the genome in the first read. Having two reads, the second is always assumed to have the opposite
# orientation of the first.
is_reverse_complement=true

# If two alignments are classified as chimeric (see min_distance), but both get assigned to the same
# annotation, this is considered a self chimera.
allow_self_chimeras=false

# Tolerance for the distance between alignments on the same read to be classified as ligated. Alignments
# from different reads have a distance of infinity.
max_ligation_distance=3


#ANNOTATION PARAMETERS

# Set the keys for the variables used for naming an annotation coming from the annotation file. Multiple
# keys can be supplied for fallback purposes in the supplied order.
name_keys=["Name", "ID"]

# Define the annotation type for sRNAs in the annotation file passed to the analysis. This relates to
# the third column in .gff files.
srna_type="sRNA"

# Since sRNA annotations can overlap with CDS or UTR annotation, they can be prioritized when annotating
# the alignments. The minimum overlap of the alignment with an annotation of type srna_type has to be at
# least min_prioritize_overlap
prioritize_srna=true
min_prioritize_overlap=0.9

# The cds_type is important if UTRs should be autocompleted or merged and for visualization purposes.
# All additional_type types are also loaded from the annotation file, if found, but not further processed.
# filter_types are also loaded without further processing and interactions with at least one partner
# belonging to one of the filter_types are excluded from the results.
cds_type="CDS"
rrna_type="rRNA"
trna_type="tRNA"
additional_types=["MySpecialType"]
filter_types=["rRNA", "tRNA"]

# If true, for every cds_type annotation, the regions upstream and downstream are set to be UTRs of types
# fiveutr_type and threeutr_type and a maximum length of autocomplete_utr_length.
autocomplete_utrs=true
autocomplete_utr_length=200
fiveutr_type="5UTR"
threeutr_type="3UTR"

# If true, all annotations of type cds_type, fiveutr_type and threeutr_type are merged to one annotation.
# All information on relative positions are with respect to the merged annotation then.
merge_utrs_and_cds=true
merge_type="CDS_UTRS"

# All space on the genome without annotation is annotated as IGR of type igr_type with the closest upstream
# and downstream annotation on the same strand used to name the IGR by separating them by a colon.
autocomplete_igrs=true
igr_type="IGR"


#FISHER EXACT TEST PARAMETERS

# For the test, a contingency table is created for every interaction. include_orientation toggles, if the
# order of the two alignments on the insert (RNA1 and RNA2) should be considered. If false, the counts of
# 2 interactions with the same partners but different orientation will be pooled together and both will
# have the same p-value assinged. include_singles decides, if non-chimeric reads assinged to an annotation
# will be included in the contingency table
include_orientation=true
include_singles=true

# Fisher's exact test uses the hypergeometric distribution to determine p-values from a contingency table
# fisher_exact_tail can be set to "left", "both" or "right". This determines, if significance is assigned
# to interactions with less occurance than to be expected by chance ("left"), more occurence than to be
# expected ("right"), or both ("both").
fisher_exact_tail="both"


#BASEPAIRING PREDICTION PARAMETERS

# Simple basepairing predictions based on alignments are performed within a region around the ligation point.
# This region is defined by the two values in bp_distance. The first one defines, how far from the ligation
# point in the direction that is "on" the read of RNA1 (upstream) and RNA2 (downstream). The second value
# corresponds to the direction that is cut off by the ligation (downstream for RNA1, upstream for RNA2) and
# should be negative to increase the distance. (45,-10) spans a region of 55 nucleotides around the ligation
# point.
bp_interval=(30,0)
bp_shift_weight=1.0

# choose, which method to use for joining pvalues from all ligation points of a given pair of annotations.
# Possible options: "stouffer", "fisher" and "fdr". "fdr" takes the minimum FDR of the pvalues.
combine_pvalues_method="stouffer"

# Set the number of sample pairs of length bp_distance_behind + bp_distance_before to be used to build the
# null model of the basepairing test.
n_genome_samples=500000

# Set the parameters for the basepairing prediction. ChimericFragments uses an alignment based approach
# that can be parameterized. The prediction is based on a substitution matrix that generates scores for
# possible alignments based on the following nucleotide pairing scores and the mismatch and gapping penalties.
AU_score=4
GC_score=5
GU_score=0
bp_mismatch_penalty=7
bp_gap_open_penalty=8
bp_gap_extend_penalty=3


#DATA FILTER PARAMETERS

# For large datasets, these parameters can be used to limit the results saved to the tables or loaded into
# the visualization. The method of Benjamini-Hochberg is applied to compute the false discovery rate.
# Set the max_fisher_fdr to filter the results by fisher's exact test. max_bp_fdr can be used to filter by the
# significance of the basepairing predictions. And interactions with a read count less than or equal to
# min_reads will be excluded from the results.
max_fisher_fdr=1.0
max_bp_fdr=1.0
min_reads=2

# toggle, if interactions without ligation points are kept in the results. Depending on the read length. With
# a sufficiently large read length, filtering out interactions without ligation points can increase the quality
# of the dataset.
keep_ints_without_ligation=true

# filter out interactions of names containing any of the specified filter queries. This option can help to
# clean mapping artifacts caused by multiple occurences of the same sequence in the genome. E.g. single reads
# from regions around rRNA or tRNA annotated as IGR can lead to many interactions in the dataset. These can be
# filtered out by specifying strings contained in the names of those annotations, e.g. ["23S", "tRNA"].
filter_name_queries = []

#DASH APP PARAMETERS

# Set the ip address and port the web server will be bound to. 0.0.0.0 will make it listen on the local
# network IP. To access the visualization, visit http://localhost:8050 in an internet browser.
address="0.0.0.0"
port=8050
