# Chimeric fragments analysis tool

This is a tool for mapping and analyzing sequencing data from RNA-Seq experiments that produce 
chimeric reads. It also provides a visualization of the results.

## Dependencies

The following procedure was tested on Ubuntu 22.04 and Arch linux (July 2022).

The script is written in the Julia programming language and requires Julia (1.7.3+) to be installed 
on your system. You can download Julia here: https://julialang.org/downloads/

The script also depends on bwa-mem2 (2.2.1+) and samtools (1.16.1+) for sequence alignments and
.bam file processing. You can find the two tools on github:

https://github.com/bwa-mem2/bwa-mem2/releases

https://github.com/samtools/samtools/releases

The default parameters of this script expect samtools and bwa-mem2 installed on the system or set
in the $PATH variable. ChimericFragments depends on a number of Julia packages that have to be 
installed and precompiled. To do this, run

>julia install.jl

The installation process could take up to a few minutes.

## Analysis procedure

The workflow can be split into the following steps:

1. create a project folder
2. add sequencing files, an annotation file and a genome sequence file
3. copy default config.jl to the folder and edit as necessary
4. run analyze.jl to analyze your experiment
5. run visualize.jl to visualize the results

### edit config.jl

IT IS IMPORTANT TO KEEP THE STRUCTURE OF config.jl INTACT. PLEASE ONLY EDIT TEXT BETWEEN "", NUMBERS
OR EXCHANGE true AND false.

All relevant analysis parameter and the dependencies can be set in the config.jl file. Some 
parameters have to be set, others are optional and affect the analysis outcome.

#### Necessary parameters:

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

#### Chimeric fragments parameters:

Here, the way how alignments are processed and chimeric fragments are classified, can be set.
min\_distance defines the minimum distance between two aligned fragments on the same strand from to
be classified as chimeric. is\_reverse\_complement has to be set, if the reads come from complementary
DNA as in the RIL-Seq experiment.

#### Annotation parameters:

The RILSeq analysis requires a complete annotation of the genome, i.e. all regions without an
annotation are ignored. The analysis can automatically generate a complete annotation, or expects
the annotation to be supplied by the annotation\_file. This can be set using the autocomple\_utrs
option in the config.jl file. If set to false, the types of the supplied UTRs and IGRs have to be
the same as defined by the options fiveutr\_type, threeutr\_type and igr\_type.

#### Fisher Exact Test parameters:

include\_read\_identity defines, if the test should take the order of the chimeric pair's placement
on the read(s) into account and include\_singles the addition of single transcripts of each partner.
max\_fdr sets the maximum false discovery rate and min\_reads defines a cut off for the output into
the final tables.

#### Visualization parameters:

The visualization of the results is implemented as a Dash web application. When it runs, it can be
accessed from within your internet browser, or any internet browser within your local network. The
parameter address and port define where the server runs. The default values will make it possible
to reach the application by visiting http://localhost:8050 in your browser.

### run analyze.jl and visualize.jl

The scripts can be executed by passing them to the julia executable. In a terminal, run

>julia path/to/analyze.jl path/to/project_folder/config.jl

To rerun the analysis (e.g. with a different set of parameters), remove the results folder from
the project folder and rerun the above command. After a successful run of the analysis, the
visualization can be started by running

>julia path/to/visualize.jl path/to/project_folder/config.jl

## Analysis results

The script will generate 2 folders called data\_processing and results. In data\_processing, 
extracted reads (fastq.gz) and alignments (.bam, .bam.bai) are saved. Also alignment stats are 
generated (.bam.log). In the results folder, several subfolder are generated and all relevant 
results are combined in two tables with filenames interactions.xlsx and singles.xlsx

## Example analysis data

To run the analysis on the supplied example data, open a terminal in this folder and run

>julia analyze.jl example_data/config.jl

The example analysis will take a minute and the resulting tables will be located in the folder 
example_data/results. Since this toy data only contains a few thousand reads, the tables should 
show between 2 and 114 interactions across the 6 conditions. 

## Copyright notice and disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. You should have received a copy of the 
GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Author: Malte Siemers, Friedrich Schiller University Jena
