#!/bin/bash
usage(){
   echo "## This script maps your reference genome to your target set using either blastn (default) or tblastx
   ## Dependencies: BLAST+
   ## ARGUMENTS
  
   ##--- The following are REQUIRED. 
   #  -R reference genome
   #  -T target fasta in nucleotide format

   ##--- The following arguments are OPTIONAL.
   #  -B what type of BLAST to use? Options currently are blastn or tblastx. The latter is potentially more sensitive but also much more time-consuming. Default = blastn.
   #  -n number of threads for blast to use. default = 1.
   #  -e evalue cutoff to retain blast matches. default = 1e-6.
   "
}
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi


BLAST_TYPE=blastn
threads=1
evalue=1e-6
while getopts R:T:n:B:e: option
do
case "${option}"
in

R) REF=${OPTARG};;
T) TARGET=${OPTARG};;
n) threads=${OPTARG};;
B) BLAST_TYPE=${OPTARG};;
e) evalue=${OPTARG};;

esac
done

REF_BASE=`basename ${REF}`
db_exist_count=`ls ${REF_BASE}_db.* | wc -l`

if [[ ${db_exist_count} == 7 ]]; then
    echo "Assuming blast database exists."
else
    echo "Making blast database..."
    makeblastdb -in ${REF} -dbtype nucl -out ${REF_BASE}_db
fi

if [[ ${BLAST_TYPE} == "tblastx" ]]; then
     echo "Using tblastx."
     tblastx -query ${TARGET} \
        -db ${REF_BASE}_db \
        -out `basename ${TARGET}`_to_${REF_BASE}.tblastx.txt.temp \
        -evalue ${evalue} \
        -num_threads ${threads} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
     echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
        cat - `basename ${TARGET}`_to_${REF_BASE}.tblastx.txt.temp > `basename ${TARGET}`_to_${REF_BASE}.tblastx.txt
     rm `basename ${TARGET}`_to_${REF_BASE}.tblastx.txt.temp
   else
     echo "Using blastn."
     blastn -query ${TARGET} \
        -db ${REF_BASE}_db \
        -out `basename ${TARGET}`_to_${REF_BASE}.blastn.txt.temp \
        -evalue ${evalue} \
        -num_threads ${threads} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"      
     echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
        cat - `basename ${TARGET}`_to_${REF_BASE}.blastn.txt.temp > `basename ${TARGET}`_to_${REF_BASE}.blastn.txt
     rm `basename ${TARGET}`_to_${REF_BASE}.blastn.txt.temp
   fi



echo "Done!"