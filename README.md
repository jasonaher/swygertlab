# swygertlab
The following is a pipeline designed for usage in data analysis of Micro-C data. 

A couple packages are needed to run this code including:
-python version 3.10+
-cooltools
-pairtools
-numpy
-matplotlib
-deeptools

The main focal point of the pipeline begins with MicC_pipeline_V3.sh and if properly working should be the only interface needed. 
This code was designed with simplicity for the user in mind, and in parallel with this running the MicC_pipeline_V3.sh with your .fastq files
should be all that is needed. 

example: sbatch MicC_pipeline_V3.sh your_microc_dataR1.fastq your_microc_dataR2.fastq

distance_decay.py is called inside the pipeline and will output a graph containing the average short contacts of the different orientations.

overlap_dist_decay.py is not called inside the pipeline and is there if you need to compare two different micro-c data sets. 
