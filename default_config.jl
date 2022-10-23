#Necessary parameters

bwa_mem2_bin = "bwa-mem2"
samtools_bin = "samtools"

annotation_file = "annotation.gff"
genome_file = "genome.fa"
data_folder = "data"

is_paired_end = true
is_single_file_paired_end = false

suffix_read1 = "_1"
suffix_read2 = "_2"

file_type = ".fastq.gz"

samplename_condition = [
    ("conditionA_rep1", "conditionA"),
    ("conditionA_rep2", "conditionA"),
    ("conditionA_rep3", "conditionA"),
    ("conditionB_rep1", "conditionB"),
    ("conditionB_rep2", "conditionB"),
    ("conditionB_rep3", "conditionB"),
    ("conditionC_rep1", "conditionC"),
    ("conditionC_rep2", "conditionC"),
    ("conditionC_rep3", "conditionC"),
]

#Chimeric fragments parameters

min_distance=1000
is_reverse_complement=true
allow_self_chimeras=false
max_ligation_distance=5
position_distribution_bins = 50

#Annotation parameters

srna_type="sRNA"
prioritize_srna=true
min_prioritize_overlap=0.9
cds_type="CDS"
additional_types=["MySpecialType"]
rrna_type="rRNA"
trna_type="tRNA"
filter_types=["rRNA", "tRNA"]
autocomplete_utrs=true
autocomplete_utr_length=200
fiveutr_type="5UTR"
threeutr_type="3UTR"
merge_utrs_and_cds=true
autocomplete_igrs=true
igr_type="IGR"
name_keys=["Name", "ID"]

#Fisher Exact Test parameters

include_read_identity=true
include_singles=true
max_fdr=0.05
min_reads=1

#Visualization parameters

address="0.0.0.0"
port=8050
