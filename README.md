# swygertlab
The following is a pipeline designed for usage in data analysis of Micro-C data. 

A couple packages are needed to run this code including:
- python version 3.10+
- cooltools
- pairtools
- numpy
- matplotlib
- deeptools

The main focal point of the pipeline begins with MicC_pipeline_V3.sh and if properly working should be the only interface needed. 
This code was designed with simplicity for the user in mind, and in parallel with this running the MicC_pipeline_V3.sh with your .fastq files
should be all that is needed. 

example: sbatch MicC_pipeline_V3.sh your_microc_dataR1.fastq your_microc_dataR2.fastq pick_a_title_for_your_data

distance_decay.py is called inside the pipeline and will output a graph containing the average short contacts of the different orientations.

overlap_dist_decay.py is not called inside the pipeline and is there if you need to compare two different micro-c data sets. 

There are a couple assumptions made by the MicC_pipeline_V3.sh pipeline mostly about the absolute paths of files, and can be changed in the 
batch script itself if desired or simply follow the file directory scheme.
Firstly, the genome of your organism is in a folder called Saccer3. # See line 30
Secondly, the size of each chromosome is in a folder called microc and the file is called sacCer3.chrSizes. # See lines 35 and 80-83

In order, for the java jar commands to work you need to download a jar file from the aiden lab github here -> https://github.com/aidenlab/juicer/wiki/Download 
this file is also assumed to be in the working directory. 

The MicC_pipeline_V3.sh pipeline also assumes that the associated python scripts are in the working directory. 
These are part of this repository and are:
- filter_pairs.py
- distance_decay.py
