#!/bin/bash

## this script maps your WGS paired-end reads to your target set and computes the coverage at each location
## Dependencies: BWA and samtools 

while getopts F:R:T:n: option
do
case "${option}"
in

F) READS1=${OPTARG};;
R) READS2=${OPTARG};;
T) TARGET=${OPTARG};;
n) threads=${OPTARG};;

esac
done

bwa index $TARGET
bwa mem -t $threads $TARGET $READS1 $READS2 | samtools view -@ $threads -F 4 -q 20 -b | samtools sort > `basename $READS1 .fastq.gz`_to_${TARGET}.sorted.bam
BAM=`basename $READS1 .fastq.gz`_to_${TARGET}.sorted.bam
samtools index $BAM

samtools depth -a $BAM > ${BAM}.coverage
