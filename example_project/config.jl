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
    ("wt_lcd_1", "wt_lcd"),
    ("wt_lcd_2", "wt_lcd"),
    ("wt_hcd_1", "wt_hcd"),
    ("wt_hcd_2", "wt_hcd"),
    ("hfq_lcd_1", "hfq_lcd"),
    ("hfq_lcd_2", "hfq_lcd"),
    ("hfq_hcd_1", "hfq_hcd"),
    ("hfq_hcd_2", "hfq_hcd"),
    ("hfq_0.2_1", "hfq_0.2"),
    ("hfq_0.2_2", "hfq_0.2"),
]

#Chimeric fragments parameter

min_distance=1000
is_reverse_complement=true
allow_self_chimeras=false

#Annotation parameter

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
autocomplete_igrs=true
igr_type="IGR"

#Fisher Exact Test parameter

include_read_identity=true
include_singles=true
max_fdr=0.05
min_reads=10
