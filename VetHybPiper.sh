#!/bin/bash
usage(){
   echo "## Bash script to use TargetVet to detect paralogs from data for many samples that have been assembled by HybPiper.
   ## It also analyses and reports patterns of missingness, and plots all of these results in various ways.
   ## ARGUMENTS
  
   ##--- The following are REQUIRED. 
   #  -V **absolute** directory of VetTargets source code: without the trailing "/"
   #  -D directory of HybPiper output (i.e. with separate folders for each sample and for genes therein): without trailing "/". Can be '.' if current directory.
   #  -O output directory suffix. output will be named TargetVet_results_<suffix> or TargetVet_results_deduped_<minidentity>_<suffix> if -d=TRUE.
   #  -T targets fasta (nucleotide): must be in directory specified in -D
   #  -S file listing sample names to process: must be in directory specified in -D
   #  -G file listing gene names to process: must be in directory specified in -D
   #  -M does the target fasta contain multiple copies per gene (TRUE or FALSE)? If TRUE, gene names in the target fasta must follow HybPiper convention, E.g. Artocarpus-gene001 and Morus-gene001 are the same gene. (Note: The gene names file (-G argument) must still have just the gene names, i.e. gene001, for example.)
   
   ##--- The following arguments are OPTIONAL.
   #  -L minimum length of blast matches to keep for analysis. Default = 150bp.
   #  -K minimum percent identity (pident) of blast matches to Keep for analysis. Default = 70.
   #  -I do IntronStats? default = TRUE
   #  -F force overwrite of DetectParalogs.R output? default = TRUE
   #  -C do per-chromosome/scaffold stats? default = TRUE
   #  -d deduplicate contigs before running analysis? default = FALSE. (Must have BBMap's dedupe.sh in current path. See https://sourceforge.net/projects/bbmap/)
   #  -m if -d=TRUE, min % identity threshold to identify duplicates? Default = 97 (to allow for removal of alleles.)
   #  -i file listing 'ingroup' samples. This is useful if you have several outgroup taxa, which often have different paralogy patterns (especially if they were used to design the target set). A separate paralog detection analysis will be conducted using only the ingroup samples.
   #  -p rooted tree in Newick format. If provided, will make an additional paralogy heatmap using this tree instead of the cluster dendrogram. All tip labels need to match those in the samples file.
   #  -P do plots for each sample with VetTargets_genome.R? Time-consuming and not necessary. Default = FALSE.
   #  -B what type of BLAST to use? Options currently are blastn or tblastx. The latter is potentially more sensitive but also much more time-consuming. Default = blastn.
   "
}
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi


## set defaults
DO_INTRON=TRUE
DO_PER_CHROM=TRUE
DEDUPE=FALSE
MINIDENTITY=97
LENGTH=150
PIDENT=70
FORCE=TRUE
DO_PLOTS=FALSE
BLAST_TYPE=blastn
PHYLO=NULL
INGROUP=NULL
## parse args
while getopts V:D:T:S:G:L:K:I:C:d:m:O:M:F:i:p:P:B: option
do
case "${option}"
in

V) VETDIR=${OPTARG};;
D) DIR=${OPTARG};;
T) TARGETS=${OPTARG};;
S) SAMPLES=${OPTARG};;
G) GENES=${OPTARG};;
L) LENGTH=${OPTARG};;
K) PIDENT=${OPTARG};;
I) DO_INTRON=${OPTARG};;
C) DO_PER_CHROM=${OPTARG};;
d) DEDUPE=${OPTARG};;
m) MINIDENTITY=${OPTARG};;
O) OUTSUFFIX=${OPTARG};;
M) MULTI=${OPTARG};;
F) FORCE=${OPTARG};;
i) INGROUP=${OPTARG};;
p) PHYLO=${OPTARG};;
P) DO_PLOTS=${OPTARG};;
B) BLAST_TYPE=${OPTARG};;

esac
done

if [[ ${DIR} == "." ]]; then
   DIR=$PWD
fi

if [[ ${DEDUPE} == "TRUE" ]]; then
   echo "Ah, so you have chosen deduplication. Summoning Bestus Bioinformaticus!\n Using MINIDENTITY=${MINIDENTITY}"
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
 echo "collating contigs for ${i}."
    while read g; do
        if [[ -f ${DIR}/${i}/${g}/${g}_contigs.fasta.gz ]]; then
            gzip -d ${DIR}/${i}/${g}/${g}_contigs.fasta.gz
        fi
        cat ${DIR}/${i}/${g}/${g}_contigs.fasta
    done < ${DIR}/${GENES} > ${CONTIGS}

    if [[ ${DEDUPE} == "TRUE" ]]; then
      dedupe.sh in=${CONTIGS} out=${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta threads=1 minidentity=${MINIDENTITY} 2> ${CONTIGS}.dedupe.log
    fi
 fi
done < ${DIR}/${SAMPLES}

# 2. blast to TARGETS
mkdir -p ${OUTDIR}/blast_out
while read i;do
 BLASTOUT=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt
 if [[ -f "$BLASTOUT" ]]; then
    echo "blast output exists for ${i}. Skipping."
 else   
      if [[ ${DEDUPE} == "TRUE" ]]; then
         SUBJECT=${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta
      else
         SUBJECT=${OUTDIR}/assemblies_collated/${i}_all_contigs.fasta
      fi

   if [[ ${BLAST_TYPE} == "tblastx" ]]; then
      echo "Using tblastx."
      tblastx -query ${TARGETS} \
         -subject ${SUBJECT} \
         -out ${BLASTOUT} \
         -evalue 1e-6 \
         -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
   else
      echo "Using blastn."
      blastn -query ${TARGETS} \
         -subject ${SUBJECT} \
         -out ${BLASTOUT} \
         -evalue 1e-6 \
         -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"      
   fi
 fi
done < ${DIR}/${SAMPLES}

# 3. run VetTargets_genome.R with --doPlots FALSE
mkdir -p ${OUTDIR}/VetTargets_genome_output
cd ${OUTDIR}/VetTargets_genome_output


while read i;do
COVSTATS=${OUTDIR}/VetTargets_genome_output/${i}_thinned_blast_result.txt
 if [[ -f "$COVSTATS" ]]; then
    echo "VetTargets_genome output exists for ${i}. Skipping."
 else
    echo "running VetTargets_genome on ${i}."
    BL=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt
    BLH=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.withHeader.txt
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
        cat - ${BL} > ${BLH}
    
    Rscript ${VETDIR}/VetTargets_genome.R --blast_file ${BLH} \
      --output_prefix ${i} \
      --min_fragment_length ${LENGTH} \
      --min_pident ${PIDENT} \
      --max_intron_length 10000 \
      --max_intron_percent 1 \
      --min_display_intron 10 \
      --doPlots ${DO_PLOTS} \
      --doIntronStats ${DO_INTRON} \
      --doCovPerChrom ${DO_PER_CHROM} \
      --multicopyTarget ${MULTI} \
      --genelist ${DIR}/${GENES} \
      --blast_type ${BLAST_TYPE}
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
   BL=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt
   cut -d'-' -f1 ${BL} | sort | uniq > targetsourcenames.txt
fi

echo "Detecting paralogs..."
if [[ ${MULTI} == "TRUE" ]]; then
   while read TSN; do
      echo "working on $TSN."
      Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} -d ${OUTDIR}/VetTargets_genome_output/${TSN} -o ${OUTDIR}/DetectParalogs_output/${TSN} -f ${FORCE} -i ${DIR}/${INGROUP} -p ${DIR}/${PHYLO}
      echo "finished $TSN."
  done < targetsourcenames.txt
else
   Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} -d ${OUTDIR}/VetTargets_genome_output -o ${OUTDIR}/DetectParalogs_output -f ${FORCE} -i ${DIR}/${INGROUP} -p ${DIR}/${PHYLO}
fi

echo "All done."