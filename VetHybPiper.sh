#!/bin/bash

# To do:
#  1. If running on cluster, automatically add local R library path to scripts. 

usage () {
   echo -e "## Bash script to use TargetVet to detect paralogs from data for many samples that have been assembled by HybPiper.
   ## It also analyses and reports patterns of missingness, and plots all of these results in various ways.
   ## ARGUMENTS
  
   ##--- The following are REQUIRED. 
   #  -V **Absolute** directory of VetTargets source code: without the trailing "/"
   #  -D Directory of HybPiper output (i.e. with separate folders for each sample and for genes therein): without the trailing "/". Can be '.' if current directory.
   #  -T Targets fasta (nucleotide)
   #  -S File listing sample names to process
   #  -G File listing gene names to process
   #  -O Output directory suffix. Output will be written to a folder in the HybPiper directory and be named TargetVet_results_<suffix> or TargetVet_results_deduped_<minidentity>_<suffix> if -d=TRUE.
   #  -M LOGICAL. Does the target fasta contain multiple copies per gene (TRUE or FALSE)? If TRUE, gene names in the target fasta must follow HybPiper convention, E.g. Artocarpus-gene001 and Morus-gene001 are the same gene. (Note: The gene names file (-G) must still have just the gene names, i.e. gene001, for example.)
   
   ##--- The following arguments are OPTIONAL.
   #  -B What type of BLAST to use? Options currently are blastn or tblastx. The latter is potentially more sensitive but also much more time-consuming. Default = blastn.
   #  -s LOGICAL. Whether to blast each target's contigs to it separately. 
   #     Much slower (use -X TRUE to speed it up) but it should be used 
         if you expect that some of your targets are very similar to each other
         (e.g. if you've been trying to assemble paralogs individually)
   #  -E What e-value to use as a cutoff for keeping matches in blast results. Default = 1e-6.
   #  -L Minimum length of blast matches to keep for analysis. Default = 150bp.
   #  -K Minimum percent identity (pident) of blast matches to Keep for analysis. Default = 70.
   #  -I LOGICAL. Do IntronStats? default = TRUE
   #  -F LOGICAL. Force overwrite of DetectParalogs.R output? default = FALSE
   #  -C LOGICAL. Do per-chromosome/scaffold stats? default = TRUE
   #  -d LOGICAL. Deduplicate contigs before running analysis? default = FALSE. (Must have BBMap's dedupe.sh in current path. See https://sourceforge.net/projects/bbmap/)
   #  -m If -d=TRUE, min % identity threshold to identify duplicates? Default = 97 (to allow for removal of alleles.)
   #  -i File listing 'ingroup' samples. This is useful if you have several outgroup taxa, which often have different paralogy patterns (especially if they were used to design the target set). A separate paralog detection analysis will be conducted using only the ingroup samples.
   #  -p Rooted phylogeny in Newick format. If provided, will make an additional paralogy heatmap using this tree instead of the cluster dendrogram. All tip labels need to match those in the samples file.
   #  -P LOGICAL. Do plots for each sample with VetTargets_genome.R? Time-consuming and not necessary. Default = FALSE.
   #  -X LOGICAL. Parallelise certain operations using GNU Parallel? Default = FALSE
   #  -t Number of threads to use for multithreaded operations? This includes GNU parallel if -X is TRUE and BLAST if -s is FALSE.
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
FORCE=FALSE
DO_PLOTS=FALSE
BLAST_TYPE=blastn
PHYLO=NULL
INGROUP=NULL
PARALLEL=FALSE
THREADS=1
EVALUE=1e-6
SEPARATE=FALSE
## parse args
while getopts V:D:T:S:G:L:K:I:C:d:m:O:M:F:i:p:P:B:X:t:E:s: option
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
E) EVALUE=${OPTARG};;
X) PARALLEL=${OPTARG};;
t) THREADS=${OPTARG};;
s) SEPARATE=${OPTARG};;

esac
done

echo "Starting VetHybPiper.sh."

# Check that things will run.

if [[ ${PARALLEL} == "TRUE" ]] && [[ $(which env_parallel.bash | wc -l) -eq 0 ]]; then
   echo "You have chosen to use parallelisation but it does not seem to be in your PATH."
   echo -e "You can add it to your path by running \nexport PATH='\$PATH:path/to/parallel'"
   echo "Exiting."
   exit
   #echo "Setting parallel to FALSE and proceeding."
   #PARALLEL=FALSE
fi

if [[ ${DEDUPE} == "TRUE" ]] && [[ $(which dedupe.sh | wc -l) -eq 0 ]]; then
   echo "You have chosen to use deduplication but dedupe.sh it does not seem to be in your PATH."
   echo -e "You can add it to your path by running \nexport PATH='\$PATH:path/to/BBMap-or-BBTools'"
   echo "Exiting."
   exit
   #echo "Setting parallel to FALSE and proceeding."
   #PARALLEL=FALSE
fi

if [[ ${PARALLEL} == "TRUE" ]] && [[ ${THREADS} -eq 1 ]]; then
   echo "You have chosen to use parallelisation but not specified more than 1 thread. Use the -t argument or specify -X FALSE. Exiting."
   exit
fi

# we need to source env_parallel.bash in order to make env_parallel a function
# this has to be done *inside* VetHybPiper.sh because it is invoked as a separate bash process to
# the one running on the compute node (assuming this is running on a cluster)
if [[ ${PARALLEL} == "TRUE" ]]; then
	. `which env_parallel.bash`
fi

# Sort out directories
if [[ ${DIR} == "." ]]; then
   DIR=$PWD
fi

if [[ ! -f ${GENES} ]] && [[ ! -f ${DIR}/${GENES} ]]; then
   echo "Couldn't find gene list ${GENES}. Exiting."
   exit
fi
if [[ ! -f ${SAMPLES} ]] && [[ ! -f ${DIR}/${SAMPLES} ]]; then
   echo "Couldn't find sample list ${SAMPLES}. Exiting."
   exit
fi
if [[ ! -f ${TARGETS} ]] && [[ ! -f ${DIR}/${TARGETS} ]]; then
   echo "Couldn't find target fasta ${TARGETS}. Exiting."
   exit
fi


# set the name of the output directory
if [[ ${DEDUPE} == "TRUE" ]]; then
  echo -e "Ah, so you have chosen deduplication. Summoning Bestus Bioinformaticus! (aka BBTools).\n Using MINIDENTITY=${MINIDENTITY}"
  OUTDIR=${DIR}/TargetVet_results_deduped_${MINIDENTITY}_${OUTSUFFIX}
else
  OUTDIR=${DIR}/TargetVet_results_${OUTSUFFIX}
fi

# (FOR NOW, DON'T UNCOMMENT THIS!) Set the name of the output directory (don't put it in the HybPiper directory)  -- This needs a lot of tweaking to make it work...
# if [[ ${DEDUPE} == "TRUE" ]]; then
   # echo -e "Ah, so you have chosen deduplication. Summoning Bestus Bioinformaticus! (aka BBTools).\n Using MINIDENTITY=${MINIDENTITY}"
   # OUTDIR=TargetVet_results_deduped_${MINIDENTITY}_${OUTSUFFIX}
# else
   # OUTDIR=TargetVet_results_${OUTSUFFIX}
# fi


echo "All output will be written to ${OUTDIR}"
mkdir -p ${OUTDIR}

# If full paths are provided and they're not pointing to files within ${DIR} subdirectories, 
# change them to basename and copy into ${DIR} if they're not already there.
is_full_path () { echo $1 | grep '/' | wc -l; }
FP=$(is_full_path ${GENES})
if [[ $FP -eq 1 ]] && [[ ! -f ${DIR}/${GENES} ]]; then
   if [[ ! -f ${DIR}/$(basename ${GENES}) ]];then
      cp ${GENES} ${DIR}
   fi
   GENES=$(basename ${GENES})
fi
FP=$(is_full_path ${SAMPLES})
if [[ $FP -eq 1 ]] && [[ ! -f ${DIR}/${SAMPLES} ]]; then
   if [[ ! -f ${DIR}/$(basename ${SAMPLES}) ]];then
      cp ${SAMPLES} ${DIR}
   fi
   SAMPLES=$(basename ${SAMPLES})
fi
FP=$(is_full_path ${TARGETS})
if [[ $FP -eq 1 ]] && [[ ! -f ${DIR}/${TARGETS} ]]; then
   if [[ ! -f ${DIR}/$(basename ${TARGETS}) ]];then
      cp ${TARGETS} ${DIR}
   fi
   TARGETS=$(basename ${TARGETS})
fi

# Fix CRLF line endings if they exist, as they mess with the while loops
CRLF=$(file ${DIR}/${GENES} | grep 'CRLF' | wc -l)
if [[ $CRLF -eq 1 ]]; then
   echo "Gene list ${DIR}/${GENES} has CRLF line endings. Creating fixed file."
   sed 's/\r$//' ${DIR}/${GENES} > ${DIR}/${GENES}.CRLF-corrected.txt 
   GENES=${GENES}.CRLF-corrected.txt 
fi
CRLF=$(file ${DIR}/${SAMPLES} | grep 'CRLF' | wc -l)
if [[ $CRLF -eq 1 ]]; then
   echo "Sample list ${DIR}/${SAMPLES} has CRLF line endings. Creating fixed file."
   sed 's/\r$//' ${DIR}/${SAMPLES} > ${DIR}/${SAMPLES}.CRLF-corrected.txt 
   SAMPLES=${SAMPLES}.CRLF-corrected.txt 
fi


if [[ ${SEPARATE} == "FALSE" ]];then
   # 1. make TargetVet results folder and collate spades contigs per sample
   mkdir -p ${OUTDIR}/assemblies_collated
   while read i;do
      CONTIGS=${OUTDIR}/assemblies_collated/${i}_all_contigs.fasta
      if [[ -f ${CONTIGS} ]]; then
         echo "Collated contigs exist for ${i}. Skipping."
      else
         echo "Collating contigs for ${i}."

         if [[ ${PARALLEL} == "TRUE" ]]; then
            ## parallel cat
            cat_contigs(){
               if [[ -f ${DIR}/${i}/${1}/${1}_contigs.fasta.gz ]]; then
                     gzip -d ${DIR}/${i}/${1}/${1}_contigs.fasta.gz
               fi
               cat ${DIR}/${i}/${1}/${1}_contigs.fasta
            }
            # export -f cat_contigs
            echo "Collating contigs in parallel for sample ${i}."
            env_parallel --env cat_contigs \
               --env DIR \
               --env i \
               -j ${THREADS} "cat_contigs {}" :::: ${DIR}/${GENES} > ${CONTIGS}
         else
            while read g; do
               if [[ -f ${DIR}/${i}/${g}/${g}_contigs.fasta.gz ]]; then
                  gzip -d ${DIR}/${i}/${g}/${g}_contigs.fasta.gz
               fi
               cat ${DIR}/${i}/${g}/${g}_contigs.fasta
            done < ${DIR}/${GENES} > ${CONTIGS}
         fi

         if [[ ${DEDUPE} == "TRUE" ]] && [[ ${PARALLEL} == "FALSE" ]]; then
            dedupe.sh in=${CONTIGS} out=${OUTDIR}/assemblies_collated/${i}_all_contigs_deduped_${MINIDENTITY}.fasta threads=1 minidentity=${MINIDENTITY} 2> ${CONTIGS}.dedupe.log
         fi
      fi
   done < ${DIR}/${SAMPLES}

   if [[ ${DEDUPE} == "TRUE" ]] && [[ ${PARALLEL} == "TRUE" ]]; then
      echo "Running dedupe.sh on collated contigs in parallel across samples."
      env_parallel \
         --env OUTDIR \
         --env MINIDENTITY \
         -j ${THREADS} "dedupe.sh in=${OUTDIR}/assemblies_collated/{}_all_contigs.fasta out=${OUTDIR}/assemblies_collated/{}_all_contigs_deduped_${MINIDENTITY}.fasta threads=1 minidentity=${MINIDENTITY} 2> ${OUTDIR}/assemblies_collated/{}_all_contigs.fasta.dedupe.log" :::: ${DIR}/${SAMPLES}
   fi

   # 2. blast to TARGETS
   mkdir -p ${OUTDIR}/blast_out
   run_blast(){
      local BLASTOUT=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${1}_all_contigs.txt
      if [[ -f ${BLASTOUT} ]]; then
         echo "BLAST output exists for ${1}. Skipping. If you wish to redo the ${BLAST_TYPE} step, move the folder 'blast_out' to another directory, rename it to something else, or remove the relevant files in the blast_out directory with the prefix ${BLAST_TYPE}."
      else   
         # if blast output doesn't exist, run blast
         if [[ ${DEDUPE} == "TRUE" ]]; then
            local SUBJECT=${OUTDIR}/assemblies_collated/${1}_all_contigs_deduped_${MINIDENTITY}.fasta
         else
            local SUBJECT=${OUTDIR}/assemblies_collated/${1}_all_contigs.fasta
         fi

         if [[ ${BLAST_TYPE} == "tblastx" ]]; then
            echo "Mapping contigs to targets using tblastx for ${1}."
            tblastx -query ${DIR}/${TARGETS} \
               -subject ${SUBJECT} \
               -out ${BLASTOUT} \
               -evalue ${EVALUE} \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
         else
            echo "Mapping contigs to targets using blastn for ${1}."
            blastn -query ${DIR}/${TARGETS} \
               -subject ${SUBJECT} \
               -out ${BLASTOUT} \
               -evalue ${EVALUE} \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"      
         fi
      fi
   }

   if [[ ${PARALLEL} == "FALSE" ]]; then
      while read i;do
         run_blast ${i}
      done < ${DIR}/${
	  }
   else
      echo "Running ${BLAST_TYPE} in parallel across samples."
      env_parallel --env DIR \
		 --env run_blast \
         --env OUTDIR \
         --env BLAST_TYPE \
         --env TARGETS \
         --env DEDUPE \
         --env MINIDENTITY \
         --env EVALUE \
         -j ${THREADS} "run_blast {}" :::: ${DIR}/${SAMPLES}
   fi
else
  ## If separate == "TRUE"
   mkdir -p ${OUTDIR}/blast_out
   echo "Running separate blast searches for each target and its contigs."
   TEMP=${OUTDIR}/temp
   mkdir -p ${TEMP}
   ## if SEPARATE=TRUE, we need to skip contig collation and run blast individually for each target
   # we need to generate a gene list of genes in the actual target file (because we might have multiple copies of each target)
   grep '^>' ${DIR}/${TARGETS} | sed 's/>//g' | sed 's/\r$//' > ${TEMP}/gene_names.txt
   echo "Writing targets to separate fastas in ${TEMP} if they don't exist."
   # pull out the targets from the target file
   # grep command from https://www.biostars.org/p/49820/
   while read T;do
      if [[ ${MULTI} == "TRUE" ]]; then
         T_ID=$(echo ${T} | cut -d'-' -f2)
      else
         T_ID=${T}
      fi 
	  if [[ ! -f ${TEMP}/${T_ID}_target.fasta ]];then
		echo ">${T}" > ${TEMP}/${T_ID}_target.fasta
		grep -A 10000 -w ${T} ${DIR}/${TARGETS} | sed -n -e '1,/>/ {/>/ !{'p''}} >> ${TEMP}/${T_ID}_target.fasta
		# awk -v seq="${T}" -v RS='>' '$1 == seq {print RS $0}' ${DIR}/${TARGETS} > ${TEMP}/${T}_target.fasta
	  else
		echo ">${T}" >> ${TEMP}/${T_ID}_target.fasta
		grep -A 10000 -w ${T} ${DIR}/${TARGETS} | sed -n -e '1,/>/ {/>/ !{'p''}} >> ${TEMP}/${T_ID}_target.fasta
	  fi
   done < ${TEMP}/gene_names.txt
   # this function needs to be called for each gene, and must
   #  1. make a temporary fasta containing only the gene's sequence in the target file
   #  2. blast the contigs assembled for that gene to the temporary fasta 
   #     and write each result to a separate temporary output file, (these will be concatenated later)
   run_blast_separately(){
      # local CURRENT_GENE_CONTIGS=${DIR}/${i}/${1}/${1}_contigs.fasta
      if [[ -f ${DIR}/${i}/${1}/${1}_contigs.fasta ]];then
         if [[ ! -f ${TEMP}/${BLAST_TYPE}_${i}_to_${1}.txt ]]; then
            if [[ ${BLAST_TYPE} == "tblastx" ]]; then
               tblastx -query ${TEMP}/${1}_target.fasta \
                  -subject ${DIR}/${i}/${1}/${1}_contigs.fasta \
                  -out ${TEMP}/${BLAST_TYPE}_${i}_to_${1}.txt \
                  -evalue ${EVALUE} \
                  -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" 
            else
               blastn -query ${TEMP}/${1}_target.fasta \
                  -subject ${DIR}/${i}/${1}/${1}_contigs.fasta \
                  -out ${TEMP}/${BLAST_TYPE}_${i}_to_${1}.txt \
                  -evalue ${EVALUE} \
                  -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" 
            fi
         fi   
      fi
   }
   # now we loop over samples and
   # 1. blast targets to their contigs separately either as a loop or in parallel, then
   # 2. concatenate the results into one blast file to pass to VetTargets_genome.R
   while read i;do
    if [[ -f ${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt ]];then
	  echo "BLAST output exists for ${i}. Skipping. If you wish to redo the ${BLAST_TYPE} step, move the folder 'blast_out' to another directory, rename it to something else, or remove the relevant files in the blast_out directory with the prefix ${BLAST_TYPE}."
    else
      echo "Mapping contigs to targets using ${BLAST_TYPE} for ${i}, separately for each target's contigs."
      if [[ ${PARALLEL} == "FALSE" ]]; then 
         while read g; do
            run_blast_separately ${g}
         done < ${DIR}/${GENES}
      else
      echo "Doing this in parallel..."
         env_parallel --env run_blast_separately \
            --env OUTDIR \
            --env BLAST_TYPE \
            --env TARGETS \
            --env TEMP \
            --env DIR \
            --env EVALUE \
            --env i \
            -j ${THREADS} "run_blast_separately {}" :::: ${DIR}/${GENES}
      fi
      # concatenate results
      cat ${TEMP}/${BLAST_TYPE}_${i}_to_*.txt \
         > ${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt
    fi
   done < ${DIR}/${SAMPLES}
# rm -rf ${TEMP}
fi

# 3. run VetTargets_genome.R
mkdir -p ${OUTDIR}/VetTargets_genome_output
cd ${OUTDIR}/VetTargets_genome_output
run_VetTargets(){
      if [[ -f ${OUTDIR}/VetTargets_genome_output/${1}_CoverageStats_AcrossChromosomes.txt ]]; then
         echo "VetTargets_genome coverage stats output exists for ${1}. Skipping."
      else
         echo "Running VetTargets_genome on ${1}."
         local BL=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${1}_all_contigs.txt
         local BLH=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${1}_all_contigs.withHeader.txt
         echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore" | \
            cat - ${BL} > ${BLH}
         
         Rscript ${VETDIR}/VetTargets_genome.R --blast_file ${BLH} \
            --output_prefix ${1} \
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
   }
if [[ ${PARALLEL} == "FALSE" ]]; then
   while read i;do
      run_VetTargets ${i}
   done < ${DIR}/${SAMPLES}
else
   echo "Running VetTargets_genome.R in parallel across samples."
   env_parallel --env run_VetTargets \
      --env OUTDIR \
      --env BLAST_TYPE \
      --env TARGETS \
      --env VETDIR \
      --env LENGTH \
      --env PIDENT \
      --env DO_PLOTS \
      --env DO_INTRON \
      --env DO_PER_CHROM \
      --env MULTI \
      --env DIR \
      --env GENES \
      -j ${THREADS} "run_VetTargets {}" :::: ${DIR}/${SAMPLES}
fi

# 4. collate CoverageStats_AcrossChromosomes.txt
# AND
# 5. detect paralogs using paralogy Rate and step-function/nplr analysis
# AND
# 6. Output genelists for paralogs and single-copy

if [[ ${MULTI} == "TRUE" ]]; then
   # get prefixes so as to check for pre-existing covstats and run DetectParalogs.R on each separately
   i=`head -n1 ${DIR}/${SAMPLES}`
   BL=${OUTDIR}/blast_out/${BLAST_TYPE}_`basename ${TARGETS}`_to_${i}_all_contigs.txt
   cut -d'-' -f1 ${BL} | sort | uniq > targetsourcenames.txt
fi

if [[ -f ${INGROUP} ]];then
   echo "Ingroup file found."
   ING="-i ${INGROUP}"
else
   echo "Ingroup file not found or not provided."
   ING=""
fi

if [[ -f ${PHYLO} ]];then
   echo "Phylogeny file found."
   PHY="-p ${PHYLO}"
else
   echo "Phylogeny file not found or not provided."
   PHY=""
fi



echo "Detecting paralogs..."
if [[ ${MULTI} == "TRUE" ]]; then
   while read TSN; do
      echo "Working on $TSN..."
      Rscript ${VETDIR}/DetectParalogs.R -s ${DIR}/${SAMPLES} \
         -d ${OUTDIR}/VetTargets_genome_output/${TSN} \
         -o ${OUTDIR}/DetectParalogs_output/${TSN} \
         -f ${FORCE} \
         ${ING} \
         ${PHY}
      echo "Finished $TSN."
  done < targetsourcenames.txt
else
   Rscript ${VETDIR}/DetectParalogs.R \
      -s ${DIR}/${SAMPLES} \
      -d ${OUTDIR}/VetTargets_genome_output \
      -o ${OUTDIR}/DetectParalogs_output \
      -f ${FORCE} \
      ${ING} \
      ${PHY}
fi

echo "VetHybPiper.sh finished!"