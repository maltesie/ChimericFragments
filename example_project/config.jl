#Necessary parameter

bwa_mem2_bin = "bwa-mem2"
samtools_bin = "samtools"

annotation_file = "annotation.gff"
genome_file = "genome.fa"
data_folder = "data"

is_paired_end = true

suffix_read1 = "_1"
suffix_read2 = "_2"

samplename_condition = [
    ("bc1", "conditionA"),
    ("bc2", "conditionA"),
    ("bc3", "conditionB"),
    ("bc4", "conditionB"),
    ("bc5", "conditionC"),
    ("bc6", "conditionC"),
]

#Chimeric fragments parameter

min_distance=1000
is_reverse_complement=true

#Annotation parameter

srna_type="sRNA"
cds_type="CDS"
additional_types = []
rrna_type="rRNA"
trna_type="tRNA"
filter_types=["rRNA", "tRNA"]
autocomplete_utrs=true
autocomplete_utr_length=200
fiveutr_type="5UTR"
threeutr_type="3UTR"
autocomplete_igrs=true
igr_type="IGR"

#Fisher Exact Test parameter

include_read_identity=true
include_singles=true
max_fdr=0.05
min_reads=1
