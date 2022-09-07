# Chimeric fragments analysis tool

This is a configurable script for demultiplexing, mapping and analyzing sequencing data from RNA-Seq 
experiments that produce chimeric reads. This document serves as an installation guide and manual.

## Dependencies

The following procedure was tested on Ubuntu 22.04 and Arch linux.

The script is written in the Julia programming language and requires Julia (1.7.3+) to be installed 
on your system. You can download Julia here: https://julialang.org/downloads/

The script also depends on bwa-mem2 and samtools for sequence alignments .bam file generation. You
can find the two tools on github:

https://github.com/bwa-mem2/bwa-mem2/releases

https://github.com/samtools/samtools/releases

The default parameters of this script expect samtools and bwa-mem2 installed on the system or set
in the $PATH variable.

## Analysis procedure

The workflow can be split into a few simple steps:

1. create a project folder
2. copy 2 paired end sequencing files, an annotation file and a genome sequence file into it
3. copy default config.jl into project folder and edit as necessary
4. run the chimeric_fragments.jl script and pass the project folder as an argument

### edit config.jl

IT IS IMPORTANT TO KEEP THE STRUCTURE OF config.jl INTACT. PLEASE ONLY EDIT TEXT BETWEEN "", NUMBERS
OR EXCHANGE true AND false.

All relevant analysis parameter and the dependencies can be set in the config.jl file. Some 
parameters have to be set, others are optional and affect the analysis outcome.

#### Necessary parameter:

If bwa-mem2 and samtools are installed systemwide, the options bwa\_mem2\_bin and samtools\_bin can be
left unchanged, otherwise they should point to the executables of bwa-mem2 and samtools.

read1\_file, read2\_file, annotation\_file, genome\_file have to be set according to the file names
of the respective file. The files have to be in the project folder next to the config.jl. The
annotation file has to be in the GFF3 format (.gff, .gff3), the reference genome sequence has to be
supplied as a fasta compatible file (.fa, .fasta, .fna). The read files can be eather fasta or fastq
compatible and can be gzipped or not.

samplename\_condition\_antibarcode defines each sample in the RILseq experiment by assigning the
corresponding barcode to a sample name and its condition. For every sample, one line in () has to be
set, providing a sample name, the condition and the antibarcode, each in quotes ("") and separated by 
comma as in the default config.jl.

#### Chimeric fragments parameter:

Here, the way how alignments are processed and chimeric fragments are classified, can be set.
min\_distance defines the minimum distance between two aligned fragments on the same strand from to
be classified as chimeric. is\_reverse\_complement has to be set, if the reads come from complementary
DNA as in the RIL-Seq experiment.

#### Annotation parameter:

The RILSeq analysis requires a complete annotation of the genome, i.e. all regions without an
annotation are ignored. The analysis can automatically generate a complete annotation, or expects
the annotation to be supplied by the annotation\_file. This can be set using the autocomple\_utrs
option in the config.jl file. If set to false, the types of the supplied UTRs and IGRs have to be
the same as defined by the options fiveutr\_type, threeutr\_type and igr\_type.

#### Fisher Exact Test parameter:

include\_read\_identity defines, if the test should take the order of the chimeric pair's placement
on the read(s) into account and include\_singles the addition of single transcripts of each partner.
max\_fdr sets the maximum false discovery rate and min\_reads defines a cut off for the output into
the final tables.

### run chimeric_fragments.jl

The script can be executed by passing it to the julia executable. In a terminal, run

>julia path/to/chimeric_fragments.jl path/to/project_folder

To rerun the analysis (e.g. with a different set of parameters), remove the results folder from
the project folder and rerun the above command. 

## Analysis results

The script will generate 2 folders called data\_processing and results. In data\_processing, 
extracted reads (fastq.gz) and alignments (.bam, .bam.bai) are saved. Also alignment stats are 
generated (.bam.log). In the results folder, several subfolder are generated and all relevant 
results are combined in two tables with filenames interactions.xlsx and singles.xlsx
