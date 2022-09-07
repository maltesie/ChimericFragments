#Necessary parameter

bwa_mem2_bin = "bwa-mem2"
samtools_bin = "samtools"

read1_file = "sample_r1.fastq.gz"
read2_file = "sample_r2.fastq.gz"

annotation_file = "annotation.gff"

genome_file = "genome.fa"

samplename_condition_antibarcode = [
    ("condition1_1", "condition1", "AATAATGT"),
    ("condition1_2", "condition1", "CAACACTT"),
    ("condition2_1", "condition2", "ATAATTCT"),
    ("condition2_2", "condition2", "GTCCATAT"),
    ("condition3_1", "condition3", "CAAGTGAT"),
    ("condition3_2", "condition3", "CGACTTGG"),
    ("condition4_1", "condition4", "GCGAGTTG"),
    ("condition4_2", "condition4", "AAGACGGG"),
    ("condition5_1", "condition5", "CTGTAGGG"),
    ("condition5_2", "condition5", "GCCGAGGG"),
    ("condition6_1", "condition6", "GGTACCGG"),
    ("condition6_2", "condition6", "GCCCTCCG"),
]

#Chimeric fragments parameter

min_distance=1000
is_reverse_complement=true

#Annotation parameter

srna_type="sRNA"
cds_type="CDS"
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
