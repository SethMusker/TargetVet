

TargetSupercontigs<-function(blast_file,genome,output_prefix,blacklist=NULL,btype="blastn",min_length=1000,max_length=10000,keepDuplicates=FALSE,keepLongSupercontigs=FALSE){
  
  # exit if any files don't exist
  if(any(!file.exists(genome,blast_file))){
    stop("Blast result and/or genome files don't exist. Check spelling and paths.\nExiting...\n")
  }else{
    if(!is.null(blacklist)){
      stopifnot("Blacklist file doesn't exist. Check spelling and path.\nExiting...\n" = file.exists(blacklist))
    }
    ######################
    ### START FUNCTION ###
    ######################
    
    if(btype=="blastn"){
      cat("Note: assuming blast file is output from blastn\n")
      btypefac<-1
    }else{
      if(btype=="tblastn"){
        cat("Note: assuming blast file is output from tblastn\n")
        btypefac<-3
      }
    }
    if(!file.exists(paste0(genome,".fai"))) {
      cat("No index file detected (searched for ",paste0(genome,".fai"),"). Indexing...\n",sep="")
      indexFa(genome)
      cat("done\n")
    }
    
    # read in blast file
    dat<-as_tibble(read.table(blast_file,header=T,stringsAsFactors=T))
    
    # check that all sseqids are in the genome by reading in the fai
    myindex<-read.table(paste0(genome,".fai"))
    if(any(! dat$sseqid %in% myindex[,1]) ){
      stop("Some or all sseqids in the blast file are not found in the genome's index. Have you got the right genome version? \n
            If yes, consider deleting the index and rerunning this script, which will automatically reindex using Rsamtools (equivalent to samtools faidx).\n
            Exiting...\n")
    }else{
      if(!is.null(blacklist)){
        cat("Removing targets in blacklist:", blacklist,"\n")
        black<-suppressMessages(scan(blacklist,what="character"))
        dat<-dat[! dat$qseqid %in% black,]
      }
      dat_Ranges<-dat %>% 
        group_by(qseqid,sseqid) %>%
        summarise(qlen_nucleotides=unique(qlen)*btypefac,
                  range_start=min(min(sstart),min(send)),
                  range_end=max(max(sstart),max(send)),
                  orientation=ifelse(min(sstart)>min(send),"-","+"),
                  new_target_length=range_end-range_start,
                  new_target_contig_length=unique(slen),
                  new_target_prop_contig_length=round(new_target_length/new_target_contig_length,4),
                  n_nucleotides_added=(range_end-range_start)-unique(qlen)*btypefac,
                  mean_pident=round(mean(pident),2),
                  percent_length_increase=round((new_target_length-qlen_nucleotides)/qlen_nucleotides,3)*100,
                  .groups = "keep") %>%
        filter(new_target_length >= min_length) %>%
        arrange(qseqid)
      
	  # Find and maybe remove duplicates 
      dups<-names(which(table(dat_Ranges$qseqid)>1))
      if(length(dups)>0){
        cat("Warning: after removing blacklisted targets and supercontigs shorter than ",min_length,"bp, ",length(dups)," targets still have multiple hits, suggesting paralogy.
        \nTheir names are written to file: ",output_prefix,"__warn_duplicated_TargetSupercontigs.txt\n",sep="")
        write.table(data.frame(dups=dups),paste0(output_prefix,"__warn_duplicated_TargetSupercontigs.txt"),row.names = F,col.names = F,quote = F)
        if(!keepDuplicates){
          cat("Duplicated loci will NOT be retained in the resulting supercontig target fasta.\n")
          dat_Ranges<-dat_Ranges[!dat_Ranges$qseqid %in% dups,]
        } else {
          cat("Duplicated loci will be retained in the resulting supercontig target fasta.\n")
        }
        chosen_regions<-paste0(dat_Ranges$sseqid,":",dat_Ranges$range_start,"-",dat_Ranges$range_end,":",dat_Ranges$orientation)
        chosen_Granges<-as(chosen_regions,"GRanges")
      }
      
      chosen_regions_sequences<-scanFa(genome,chosen_Granges)
      
      ## plot histogram of sequence lengths
      data.frame(Length=width(chosen_regions_sequences)) %>%
        ggplot(aes(x=Length))+
        geom_histogram(fill="skyblue",colour="black",binwidth = 500)+
        geom_vline(xintercept = max_length,colour="red",size=1.5)+
        theme_bw()+
        # scale_x_continuous(breaks = round(seq(min(width(chosen_regions_sequences)), max(width(chosen_regions_sequences)), by = 1000),1))+
        scale_x_continuous(breaks = scales::breaks_width(1000))+
        scale_y_continuous(breaks = scales::breaks_width(2))
      suppressMessages(ggsave(paste0(output_prefix,"_TargetSupercontigs_length_histogram.pdf")))
      
      ## NOTE: rename sequences with target and genome IDs
      genomeName<-basename(genome)
      New_names<-dat_Ranges %>% transmute(ID=paste0("Supercontig_",genomeName,"_",sseqid,"_region-",range_start,"-",range_end,"_Target_",qseqid))
      if(all(names(chosen_regions_sequences)==New_names$sseqid)){
        names(chosen_regions_sequences)<-New_names$ID
        # filter out long supercontigs
        if(!keepLongSupercontigs){
          cat("Removing",sum(width(chosen_regions_sequences)>=max_length),"supercontigs >= max_length:",max_length,"\n")
          chosen_regions_sequences<-chosen_regions_sequences[-which(width(chosen_regions_sequences)>=max_length)]
        }else{
          cat("Note: Ignoring max_length and keeping long supercontigs\n")
        }
        cat("Writing",length(chosen_regions_sequences),"supercontigs to", paste0(output_prefix,"_TargetSuperContigs.fasta\n"))
        writeXStringSet(chosen_regions_sequences,paste0(output_prefix,"_TargetSuperContigs.fasta"))
        # cat(width(chosen_regions_sequences),"\n")
      }else{
        stop("There was a problem appending the target names to the supercontig names. Exiting\n")
      }
    }
  }
}

#### DONE DEFINING FUNCTIONS
suppressMessages(suppressWarnings(require(optparse,quietly=TRUE,warn.conflicts=FALSE)))

p <- OptionParser(usage=" This script will take\n
1. a tabular blast result (-outfmt 6) ***WITH A HEADER LINE** (see the readme) with your target exons as the query and (draft) genome as the subject\n
2. your genome fasta\n
3. Optionally a blacklist of targets to ignore\n
  and do the following:\n
   a. Find the start and end coordinates of the supercontigs according to the blast result.\n
   b. Flag and (by defualt) remove targets with more than one inferred supercontig match. You can change this by setting keepDuplicates to TRUE.\n
   c. Remove any shorter than min_length\n
   d. Index the reference genome (if no index is present) and extract the supercontigs\n
   e. If keepLongSupercontigs=FALSE (the default), remove any supercontigs longer than max_length\n
   f. Write a Fasta with your new supercontigs\n
  Run using Rscript, e.g.\n
  Rscript TargetSupercontigs.R --blast_file blastn_targets_to_genome.txt --genome my_genome.fa --output_prefix test --blacklist targets.remove.txt\n"
)
# Add a positional argument
p <- add_option(p, c("-b","--blast_file"), help="<Required: tab-delimited blast result, target=query, genome=subject. **With Header!**>",type="character")
p <- add_option(p, c("-g","--genome"), help="<Required: genome fasta used in blast search>",type="character")
p <- add_option(p, c("-o","--output_prefix"), help="<Required: prefix to name results files>",type="character")
p <- add_option(p, c("-x","--blacklist"), help="<Optional: file listing targets to exclude>",type="character")
p <- add_option(p, c("--btype"), help="<Optional: what type of blast search? Options are blastn (default) and tblastn>",type="character",default="blastn")
p <- add_option(p, c("--min_length"), help="<Optional: supercontigs shorter are ignored; default 1kb>",type="numeric",default=1000)
p <- add_option(p, c("--max_length"), help="<Optional: supercontigs longer are ignored; default 10kb>",type="numeric",default=10000)
p <- add_option(p, c("--keepDuplicates"), help="<Optional: logical: whether to keep targets with more than one hit in the genome; default FALSE>",type="logical",default=FALSE)
p <- add_option(p, c("--keepLongSupercontigs"), help="<Optional: logical: whether to keep supoercontigs longer than max_length; default FALSE>",type="logical",default=FALSE)

# parse
args<-parse_args(p)

if(any(is.null(c(args$blast_file,args$genome,args$output_prefix)))) {
  # Print the help message
  print_help(p)
  cat("Hey hey hey! Some required arguments are missing.\n")
}else{
# suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dplyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(Biostrings,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(Rsamtools,quietly=TRUE,warn.conflicts=FALSE)))
## RUN
TargetSupercontigs(blast_file=args$blast_file,
                   genome=args$genome,
                   output_prefix=args$output_prefix,
                   blacklist=args$blacklist,
                   btype=args$btype,
                   min_length=args$min_length,
                   max_length=args$max_length,
                   keepDuplicates=args$keepDuplicates,
                   keepLongSupercontigs=args$keepLongSupercontigs)
}