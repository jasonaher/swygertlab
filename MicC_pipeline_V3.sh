#!/bin/bash 

#SBATCH --partition=amilan
#SBATCH --job-name=micc_pipeline
#SBATCH --output=test_micro-c_pipeline.%j.out
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --nodes=2
#SBATCH --ntasks=18
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasonher@colostate.edu

module purge
module load anaconda
conda activate MicC # you will need an active conda environment with the bowtie2, samtools, pairtools, and cooler packages 

r1=$1 
r2=$2
#file_names=$3 # include how you want the files to be named 

sample1=$(echo ${r1} | sed 's/\..*$//')
sample2=$(echo ${r2} | sed 's/\..*$//')
#sample3=${file_names}
sample3=$sample1$sample2
#bzip2 -d ${r1}
#bzip2 -d ${r2}

#bowtie2 will align the sample reads to the reference genome
#after '-x' should be the directory that leads to files with your organism's whole genome
bowtie2 --very-sensitive -p 16 --reorder -x /pl/active/swygertlab/jasonher/Saccer3/Saccer3 -1 ${sample1}.fastq -2 ${sample2}.fastq -S ${sample3}.sam

#turning the bowtie2 .sam file output into a .bam file and using the .bam file to make a .pairs file
#in the 2nd command here after '-c' should be a .chrSizes file which contains the base pair size of each chromosome
samtools view -S -b ${sample3}.sam > ${sample3}.bam
samtools view -h ${sample3}.bam | pairtools parse -c /pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes -o ${sample3}_parsed.pairs.gz

#pairtools sort puts the reads in base pair sequential order
#dedup removes duplicates
#select removes based on some criteria here we only want unrescued reads for more info check out the pairtools select documentation
pairtools sort --nproc 8 --tmpdir=./ -o ${sample3}_sorted.pairs.gz ${sample3}_parsed.pairs.gz
pairtools dedup --mark-dups -o ${sample3}_deduped.pairs.gz ${sample3}_sorted.pairs.gz
pairtools select '(pair_type == "UU")' -o ${sample3}_filtered.pairs.gz ${sample3}_deduped.pairs.gz 
pairtools split --output-pairs ${sample3}_output.pairs.gz ${sample3}_filtered.pairs.gz

#in order to access files in python script unzipping them
#the filter_orientation_heading.py python script will generate IN, OUT, SAME, and NoFilter .pairs files the python scripts must be in working directory with command as written
gunzip ${sample3}_output.pairs.gz
python filter_orientations_heading.py ${sample3}_output.pairs
in=$(wc -l ${sample3}_output_IN_reads.pairs)
out=$(wc -l ${sample3}_output_OUT_reads.pairs)
same=$(wc -l ${sample3}_output_SAME_reads.pairs)
noIN=$(wc -l ${sample3}_output_noIN.pairs)
in_reads=$(echo ${in} | cut -d ' ' -f 1)
out_reads=$(echo ${out} | cut -d ' ' -f 1)
same_reads=$(echo ${same} | cut -d ' ' -f 1)
noIN_reads=$(echo ${noIN} | cut -d ' ' -f 1)
sum=$((${in_reads}+${out_reads}+${same_reads}))
#python distance_decay.py script generates short distance decay plots from 0 to 2000 bp
python distance_decay.py ${sample3}_output_IN_reads.pairs ${in_reads} ${sum}
python distance_decay.py ${sample3}_output_OUT_reads.pairs ${out_reads} ${sum}
python distance_decay.py ${sample3}_output_SAME_reads.pairs ${same_reads} ${sum}
python distance_decay.py ${sample3}_output_noIN.pairs ${noIN_reads} ${sum}

#rezipping files
bgzip ${sample3}_output_IN_reads.pairs
bgzip ${sample3}_output_OUT_reads.pairs
bgzip ${sample3}_output_SAME_reads.pairs
bgzip ${sample3}_output_noIN.pairs

#creating cooler files
pairix -f ${sample3}_output_IN_reads.pairs.gz
pairix -f ${sample3}_output_OUT_reads.pairs.gz
pairix -f ${sample3}_output_SAME_reads.pairs.gz
pairix -f ${sample3}_output_noIN.pairs.gz

cooler cload pairix /pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes:150 ${sample3}_output_IN_reads.pairs.gz ${sample3}_output_IN_reads.cool
cooler cload pairix /pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes:150 ${sample3}_output_OUT_reads.pairs.gz ${sample3}_output_OUT_reads.cool
cooler cload pairix /pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes:150 ${sample3}_output_SAME_reads.pairs.gz ${sample3}_output_SAME_reads.cool
cooler cload pairix /pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes:150 ${sample3}_output_noIN.pairs.gz ${sample3}_output_noIN.cool

cooler balance ${sample3}_output_IN_reads.cool
cooler zoomify ${sample3}_output_IN_reads.cool
cooler balance ${sample3}_output_OUT_reads.cool
cooler zoomify ${sample3}_output_OUT_reads.cool
cooler balance ${sample3}_output_SAME_reads.cool
cooler zoomify ${sample3}_output_SAME_reads.cool
cooler balance ${sample3}_output_noIN.cool
cooler zoomify ${sample3}_output_noIN.cool

gzip $sample1.fastq
gzip $sample2.fastq
gzip $sample3.sam
