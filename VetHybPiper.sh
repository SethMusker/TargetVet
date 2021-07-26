## Bash script to use TargetVet to detect paralogs from data for many samples that have been assembled by HybPiper


## set defaults
threads=1
LENGTH=100
## parse args
while getopts V:D:T:S:G:N:L: option
do
case "${option}"
in

V) VETDIR=${OPTARG};;
D) DIR=${OPTARG};;
T) TARGETS=${OPTARG};;
S) SAMPLES=${OPTARG};;
G) GENES=${OPTARG};;
N) threads=${OPTARG};;
L) LENGTH=${OPTARG};;

esac
done

# 1. make TargetVet results folder and collate spades scaffolds per sample
mkdir -p ${DIR}/TargetVet_results/assemblies_collated
while read i;do
    while read g; do
        if [[ -f ${DIR}/${i}/${g}/${g}_contigs.fasta.gz ]]; then
            gzip -d ${DIR}/${i}/${g}/${g}_contigs.fasta.gz
        fi
        cat ${DIR}/${i}/${g}/${g}_contigs.fasta
    done < ${GENES} > ${DIR}/TargetVet_results/assemblies_collated/${i}_all_contigs.fasta
done < ${SAMPLES}

# 2. blast to TARGETS
mkdir ${DIR}/TargetVet_results/blast_out
while read i;do
    blastn -query ${TARGETS} \
        -subject ${DIR}/TargetVet_results/assemblies_collated/${i}_all_contigs.fasta \
        -out ${DIR}/TargetVet_results/blast_out/blastn_${TARGETS}_to_${i}_all_contigs.txt \
        -evalue 1e-6 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" \
        -num_threads ${threads}
done < ${SAMPLES}

# 3. run VetTargets_genome.R with --doPlots FALSE
mkdir ${DIR}/TargetVet_results/VetTargets_genome_output
cd ${DIR}/TargetVet_results/VetTargets_genome_output
while read i;do
    BL=${DIR}/TargetVet_results/blast_out/blastn_${TARGETS}_to_${i}_all_contigs.txt
    BLH=${DIR}/TargetVet_results/blast_out/blastn_${TARGETS}_to_${i}_all_contigs.withHeader.txt
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
        cat - ${BL} > ${BLH}
    
    Rscript ${VETDIR}/VetTargets_genome.R --blast_file ${BLH} \
    --output_prefix ${i} \
    --min_fragment_length ${LENGTH} \
    --min_pident 60 \
    --max_intron_length 10000 \
    --max_intron_percent 1 \
    --min_display_intron 1 \
    --doPlots FALSE > ${i}.Rout
done < ${SAMPLES} 

# 4. collate ContigStats_acrossChromosomes.txt
# AND
# 5. detect paralogs using paralogyRate and breakpoint analysis
# AND
# 6. Output genelists for paralogs and single-copy

Rscript ${VETDIR}/DetectParalogs.R -s ${SAMPLES} -d ${DIR}/TargetVet_results/VetTargets_genome_output