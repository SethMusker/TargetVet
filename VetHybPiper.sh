#!/bin/bash
usage(){
   echo "## Bash script to use TargetVet to detect paralogs from data for many samples that have been assembled by HybPiper.
## ARGUMENTS
#  -V **absolute** directory of VetTargets code: without trailing "/"
#  -D directory of HybPiper output (i.e. with separate folders for each sample and for genes therein): without trailing "/". Can be '.' if current directory.
#  -O Output directory suffix. output will be named TargetVet_results_<suffix> or TargetVet_results_deduped_<minidentity>_<suffix> if -d=TRUE.
#  -T targets fasta (nucleotide): must be in directory specified in -D
#  -S file listing sample names to process: must be in directory specified in -D
#  -G file listing gene names to process: must be in directory specified in -D
#  -L minimum length of blast matches to keep for analysis. Default 150bp.
#  -I do IntronStats? default = TRUE
#  -M does target file contain multiple copies per gene (TRUE or FALSE)? If TRUE, gene names must follow HybPiper convention, E.g. Artocarpus-gene001 and Morus-gene001 are the same gene. (Note: The gene names file (-G argument) must still have just the gene names, i.e. e.g. gene001.)
#  -F force overwrite of DetectParalogs.R output? default = TRUE
#  -C do per-chromosome/scaffold stats? default = TRUE
#  -d deduplicate contigs before running analysis? default = FALSE. (Must have bbmap's dedupe.sh in current path. See https://sourceforge.net/projects/bbmap/)
#  -m if -d=TRUE, min % identity threshold to identify duplicates? Default = 97 (to allow for removal of alleles.)
"
}
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi


## set defaults
LENGTH=100
DO_INTRON=TRUE
DO_PER_CHROM=TRUE
DEDUPE=FALSE
MINIDENTITY=97
LENGTH=150
FORCE=TRUE
## parse args
while getopts V:D:T:S:G:L:I:C:d:m:O:M:F: option
do
case "${option}"
in

V) VETDIR=${OPTARG};;
D) DIR=${OPTARG};;
T) TARGETS=${OPTARG};;
S) SAMPLES=${OPTARG};;
G) GENES=${OPTARG};;
L) LENGTH=${OPTARG};;
I) DO_INTRON=${OPTARG};;
C) DO_PER_CHROM=${OPTARG};;
d) DEDUPE=${OPTARG};;
m) MINIDENTITY=${OPTARG};;
O) OUTSUFFIX=${OPTARG};;
M) MULTI=${OPTARG};;
F) FORCE=${OPTARG};;

esac
done

if [[ ${DIR} == "." ]]; then
   DIR=$PWD
fi

if [[ ${DEDUPE} == "TRUE" ]]; then
   OUTDIR=${DIR}/TargetVet_results_deduped_${MINIDENTITY}_${OUTSUFFIX}
else
   OUTDIR=${DIR}/TargetVet_results_${OUTSUFFIX}
fi

# 1. make TargetVet results folder and collate spades scaffolds per sample
mkdir -p ${OUTDIR}/assemblies_collated
while read i;do
CONTIGS=${OUTDIR}/assemblies_collated/${i}_all_contigs.fasta
 if [[ -f "$CONTIGS" ]]; then
    echo "collated contigs exist for ${i}. Skipping."
 else
 echo "collating contigs for ${i}"
    while read g; do
        if [[ -f ${DIR}/${i}/${g}/${g}_contigs.fasta.gz ]]; then
            gzip -d ${DIR}/${i}/${g}/${g}_contigs.fasta.gz
        fi
        cat ${DIR}/${i}/${g}/${g}_contigs.fasta
    done < ${DIR}/${GENES} > ${CONTIGS}

    if [[ ${DEDUPE} == "TRUE" ]]; then
    echo "Done collating contigs for ${i}. Now deduplicating."
      dedupe.sh in=${CONTIGS} out=${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta threads=1 minidentity=${MINIDENTITY} 2> ${CONTIGS}.dedupe.log
      # rm ${CONTIGS}
      # mv ${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta ${CONTIGS}
    fi
 fi
done < ${DIR}/${SAMPLES}

# 2. blast to TARGETS
mkdir -p ${OUTDIR}/blast_out
while read i;do
 BLASTOUT=${OUTDIR}/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt
 if [[ -f "$BLASTOUT" ]]; then
    echo "blast output exists for ${i}. Skipping."
 else   
      if [[ ${DEDUPE} == "TRUE" ]]; then
         SUBJECT=${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta
      else
         SUBJECT=${OUTDIR}/assemblies_collated/${i}_all_contigs.fasta
      fi

    blastn -query ${TARGETS} \
        -subject ${SUBJECT} \
        -out ${BLASTOUT} \
        -evalue 1e-6 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
 fi
done < ${DIR}/${SAMPLES}

# 3. run VetTargets_genome.R with --doPlots FALSE
mkdir -p ${OUTDIR}/VetTargets_genome_output
cd ${OUTDIR}/VetTargets_genome_output


while read i;do
 if [[ -f "$COVSTATS" ]]; then
    echo "VetTargets_genome output exists for ${i}. Skipping."
 else
    echo "running VetTargets_genome on ${i}"
    BL=${OUTDIR}/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt
    BLH=${OUTDIR}/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.withHeader.txt
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
        cat - ${BL} > ${BLH}
    
    Rscript ${VETDIR}/VetTargets_genome.R --blast_file ${BLH} \
    --output_prefix ${i} \
    --min_fragment_length ${LENGTH} \
    --min_pident 60 \
    --max_intron_length 10000 \
    --max_intron_percent 1 \
    --min_display_intron 1 \
    --doPlots FALSE \
    --doIntronStats ${DO_INTRON} \
    --doCovPerChrom ${DO_PER_CHROM} \
    --multicopyTarget ${MULTI} \
    --genelist ${DIR}/${GENES}
 fi
done < ${DIR}/${SAMPLES}

# 4. collate ContigStats_acrossChromosomes.txt
# AND
# 5. detect paralogs using paralogyRate and breakpoint analysis
# AND
# 6. Output genelists for paralogs and single-copy

if [[ ${MULTI} == "TRUE" ]]; then
   # get prefixes so as to check for pre-existing covstats and run DetectParalogs.R on each separately
   i=`head -n1 ${DIR}/${SAMPLES}`
   BL=${OUTDIR}/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt
   cut -d'-' -f1 ${BL} | sort | uniq > targetsourcenames.txt
fi

echo "Detecting paralogs..."
if [[ ${MULTI} == "TRUE" ]]; then
   while read TSN; do
      echo "working on $TSN"
      Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} -d ${OUTDIR}/VetTargets_genome_output/${TSN} -o ${OUTDIR}/VetTargets_genome_output/${TSN}/DetectParalogs_output -f ${FORCE}
      echo "finished $TSN"
  done < targetsourcenames.txt
else
   Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} -d ${OUTDIR}/VetTargets_genome_output -o ${OUTDIR}/VetTargets_genome_output/DetectParalogs_output -f ${FORCE}
fi
echo "All done."