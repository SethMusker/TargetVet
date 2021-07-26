## Bash script to use TargetVet to detect paralogs from data for many samples that have been assembled by HybPiper
## ARGUMENTS
#  -V **absolute** directory of VetTargets code: without trailing "/"
#  -D **absolute** directory of HybPiper output (i.e. with separate folders for each sample and for genes therein): without trailing "/"
#  -T targets fasta (nucleotide): must be in directory specified in -D
#  -S file listing sample names to process: must be in directory specified in -D
#  -G file listing gene names to process: must be in directory specified in -D
#  -L minimum length of blast matches to keep for analysis
#  -I do IntronStats? default = TRUE


## set defaults
LENGTH=100
DO_INTRON=TRUE
DO_PER_CHROM=TRUE
## parse args
while getopts V:D:T:S:G:L:I:C: option
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

esac
done

# 1. make TargetVet results folder and collate spades scaffolds per sample
mkdir -p ${DIR}/TargetVet_results/assemblies_collated
while read i;do
 FILE=${DIR}/TargetVet_results/assemblies_collated/${i}_all_contigs.fasta
 if [[ -f "$FILE" ]]; then
    echo "collated contigs exist for ${i}. Skipping."
 else
    while read g; do
        if [[ -f ${DIR}/${i}/${g}/${g}_contigs.fasta.gz ]]; then
            gzip -d ${DIR}/${i}/${g}/${g}_contigs.fasta.gz
        fi
        cat ${DIR}/${i}/${g}/${g}_contigs.fasta
    done < ${DIR}/${GENES} > ${DIR}/TargetVet_results/assemblies_collated/${i}_all_contigs.fasta
 fi
done < ${DIR}/${SAMPLES}

# 2. blast to TARGETS
mkdir -p ${DIR}/TargetVet_results/blast_out
while read i;do
 FILE=${DIR}/TargetVet_results/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt
 if [[ -f "$FILE" ]]; then
    echo "blast output exists for ${i}. Skipping."
 else
    blastn -query ${TARGETS} \
        -subject ${DIR}/TargetVet_results/assemblies_collated/${i}_all_contigs.fasta \
        -out ${DIR}/TargetVet_results/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt \
        -evalue 1e-6 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
 fi
done < ${DIR}/${SAMPLES}

# 3. run VetTargets_genome.R with --doPlots FALSE
mkdir -p ${DIR}/TargetVet_results/VetTargets_genome_output
cd ${DIR}/TargetVet_results/VetTargets_genome_output
while read i;do
 FILE=${i}_CoverageStats_AcrossChromosomes.txt
 if [[ -f "$FILE" ]]; then
    echo "VetTargets_genome output exists for ${i}. Skipping."
 else
    echo "running VetTargets_genome on ${i}"
    BL=${DIR}/TargetVet_results/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.txt
    BLH=${DIR}/TargetVet_results/blast_out/blastn_`basename ${TARGETS}`_to_${i}_all_contigs.withHeader.txt
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
    --doCovPerChrom ${DO_PER_CHROM}
 fi
done < ${DIR}/${SAMPLES}

# 4. collate ContigStats_acrossChromosomes.txt
# AND
# 5. detect paralogs using paralogyRate and breakpoint analysis
# AND
# 6. Output genelists for paralogs and single-copy

echo "Detecting paralogs..."
Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} -d ${DIR}/TargetVet_results/VetTargets_genome_output
echo "All done."