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
3. copy default config.jl to the folder and edit as necessary
4. run analyze.jl on the config file to analyze your experiment
5. run visualize.jl on the config file to visualize the results

### edit config.jl

IT IS IMPORTANT TO KEEP THE STRUCTURE OF config.jl INTACT. PLEASE ONLY EDIT TEXT BETWEEN "", NUMBERS
OR EXCHANGE true AND false.

All relevant analysis parameters and the paths to the binaries can be set in the config.jl file. Some 
parameters have to be set, others are optional and affect the analysis outcome. You can find detailed
descriptions of every parameter in the config file.

### run analyze.jl

The scripts can be executed by passing them to the julia executable. In a terminal, run

>julia path/to/analyze.jl path/to/project_folder/config.jl

To rerun the analysis (e.g. with a different set of parameters), remove the results folder from
the project folder and rerun the above command. 

In the specified data folder, alignments (.bam, .bam.bai) and alignment stats (.bam.log) are saved. 
In the results folder, several subfolder are generated and all relevant results are combined in two 
tables with filenames interactions.xlsx and singles.xlsx

### run visualize.jl

After a successful run of the analysis, the visualization can be started by running

>julia path/to/visualize.jl path/to/project_folder/config.jl

## Example analysis data

To run the analysis on the supplied example toy data, open a terminal in this folder and run

>julia analyze.jl example_data/config.jl

The example analysis will take a minute and the resulting tables will be located in the folder 
example_data/results. 

## Copyright notice and disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. You should have received a copy of the 
GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Author: Malte Siemers, Friedrich Schiller University Jena
