# ChimericFragments

ChimericFragments is a tool for mapping and analyzing sequencing data from RNA-Seq experiments that
produce chimeric reads. It is split into an analysis part and a Dash-based interactive visualization.
ChimericFragments is highly configurable to enable analyses of different experimental workflows in
different organisms and the visualization aims to provide an additional layer of accessibility to the
data generated by those analyses.

## Dependencies

The following procedure was tested on Ubuntu 22.04 and Arch linux.

ChimericFragments is written in the Julia programming language and requires Julia (1.9.1) to
be installed on your system. You can download Julia [here](https://julialang.org/downloads/). It is
highly recommended to install Julia as described [here](https://julialang.org/downloads/platform/).
Using a julia installation from your OS repository or a package manager like conda can lead to
problems during the installation of ChimericFragments.

ChimericFragments depends on fastp (0.23.3+) for reads trimming and on bwa-mem2 (2.2.1+) and samtools (1.14.0+)
for sequence alignments and .bam file processing. If you run ChimericFragments on Linux, the binaries are provided
and do not have to be installed separately. For other operating systems (not tested), you can find the tools on github:

https://github.com/bwa-mem2/bwa-mem2/releases

https://github.com/samtools/samtools/releases

https://github.com/OpenGene/fastp/releases

## Install

The default parameters of ChimericFragments are set for running on Linux and no paths to the binary dependencies
have to be provided. For other operating systems, the paths to the binaries have to be set in the configuration
file. ChimericFragments depends on a number of Julia packages that have to be installed and precompiled. This is
done automatically on the first run of ChimericFragments and can take up to a few minutes.

#### Update

To update ChimericFragments, download the latest version and use the new files. In general, the analysis and visualization
are only compatible within one version and analysis results should be recomputed after an update.

## Using ChimericFragments

The workflow can be split into the following steps:

1. create a project folder
3. copy default_config.jl to the folder and edit as necessary
4. run analyze.jl on the config file to analyze your experiment
5. run visualize.jl on the config file to visualize the results

#### Edit config.jl

IT IS IMPORTANT TO KEEP THE STRUCTURE OF config.jl INTACT. PLEASE ONLY EDIT TEXT BETWEEN "",
NUMBERS OR EXCHANGE true AND false.

All relevant analysis parameters and the paths to the binaries can be set in the config.jl file.
For ChimericFragments to run properly, you have to set the FILE PARAMETERS and the ANNOTATION
PARAMETERS. You can find detailed descriptions of every parameter in the config file. A working
example project can be found in the example_project folder.

#### Run analyze.jl

The scripts can be executed by passing them to the julia executable. In a terminal, run

>julia path/to/analyze.jl path/to/project_folder/config.jl

To rerun the analysis (e.g. with a different set of parameters), remove the results folder from
the project folder and rerun the above command.

In the specified data folder, alignments (.bam, .bam.bai) and alignment stats (.bam.log) are saved.
In the results folder, several subfolder are generated and all relevant results are combined in two
tables with filenames interactions.xlsx and genes.xlsx

#### Run visualize.jl

After a successful run of the analysis, the visualization can be started by running

>julia path/to/visualize.jl path/to/project_folder/config.jl

#### Example analysis data

To run the analysis on the supplied example toy data, open a terminal in this folder and run

>julia analyze.jl example_project/config.jl

The example analysis will take a minute and the resulting tables and plots will be located in the
folder example_project/results. After the successful analysis part, you can start the visualization:

>julia visualyze.jl example_project/config.jl

## Copyright notice and disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. You should have received a copy of the
GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Author: Malte Siemers, Friedrich Schiller University Jena
