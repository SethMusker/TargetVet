#!/bin/bash

## this script maps your WGS paired-end reads to your target set and computes the coverage at each location
## Dependencies: BWA, samtools 

while getopts F:R:T:n:q: option
do
case "${option}"
in

F) READS1=${OPTARG};;
R) READS2=${OPTARG};;
T) TARGET=${OPTARG};;
n) threads=${OPTARG};;
q) map_qual=${OPTARG};;

esac
done

bwa index ${TARGET}
BAM=`basename ${READS1} .gz`_to_`basename ${TARGET}`.sorted.bam
bwa mem -t ${threads} ${TARGET} ${READS1} ${READS2} | samtools view -@ ${threads} -F 4 -q ${map_qual} -b | samtools sort > ${BAM}
samtools index ${BAM}

samtools depth -a ${BAM} > ${BAM}.coverage
