### This script will take a tabular blast result (-outfmt 6) with your target exons as the query and (draft) genome as the subject
###  and do the following:
###   1. Find potential paralogs by identifying genes with good matches from more than one contig/chromosome segments, or multiple positions within them, in the reference genome
###   2. Find potentially missing genes by identifying genes with no or few good matches in the reference genome
###   3. Find genes spanning huge introns
###  then
###   4. Create sets of simple segment plots to allow for visual inspection of blast alignments of genes failing and passing the above checks
###   5. Output files listing names of targets passing or failing each check, and a file of 'clean' targets passing all checks

### Arguments:
###  --blast_file <tab-delimited blast result, target=query, genome=subject>
###  --min_pident < % identity below which to discard blast matches prior to running checks >
###  --max_intron_length < maximum allowable length of any given intron in a gene >
###  --min_fragment_length < minimum length of a blast hit, ignore anything smaller >
###  --output_prefix < prefix to name results files >
## not yet implemented
###  --max_intron_percent < maximum allowable percentage of a supercontig's length consisting of introns >


## TO DO
# 1. Properly treat tblastx result - DONE (ThinBlastResult)
# 2. differentiate containment and redundancy, and add length-prioritisation if bitscore and pident are equivalent

round_perc<-function(x) round(x,3)*100     # round and convert to percentage

GetCoverageStats<-function(data,doCovPerChrom=TRUE){
  coverage_summary_chromosome_aware<-data.frame()
  coverage_summary_chromosome_UNaware<-data.frame()
  for(i in unique(data$qseqid)){
    temp<-data %>% 
      filter(qseqid==i)%>%
      arrange(qstart) 
    
    # per-target+chromosome combo 
    if(doCovPerChrom){
      for(m in unique(temp$sseqid)){
        temp2<-temp %>% filter(sseqid==m)
        target_seq<-seq(1:unique(temp2$qlen))
        for(f in 1:nrow(temp2)){
          target_seq<-c(target_seq,seq(temp2[f,]$qstart,temp2[f,]$qend))
        }
        mycounts<-table(target_seq)-1 #subtract 1 because we've appended to the original sequence
        paralogy_index<-mean(mycounts)
        paralogy_index_ignoreMissing<-mean(mycounts[mycounts>0])
        paralog_percent<-length(mycounts[mycounts>1])/unique(temp2$qlen)
        full_percent<-length(mycounts[mycounts>0])/unique(temp2$qlen)
        paralog_percent_ignoreMissing<-paralog_percent/full_percent
        unique_percent<-length(mycounts[mycounts==1])/unique(temp2$qlen)
        missing_percent<-length(mycounts[mycounts==0])/unique(temp2$qlen)
        coverage_summary_chromosome_aware<-rbind(coverage_summary_chromosome_aware,
                                                 data.frame(qseqid=unique(temp2$qseqid),
                                                            sseqid=unique(temp2$sseqid),
                                                            paralogy_index=round(paralogy_index,2),
                                                            paralogy_index_ignoreMissing=round(paralogy_index_ignoreMissing,2),
                                                            paralog_percent=round_perc(paralog_percent),
                                                            paralog_percent_ignoreMissing=round_perc(paralog_percent_ignoreMissing),
                                                            full_percent=round_perc(full_percent),
                                                            unique_percent=round_perc(unique_percent),
                                                            missing_percent=round_perc(missing_percent)))
      }
    }
    
    
    # per-target only (ignoring which chromosome(s) it matches to)
    qlenuniq<-unique(temp$qlen)
    target_seq2<-seq(1:qlenuniq)
    for(f in 1:nrow(temp)){
      target_seq2<-c(target_seq2,seq(temp[f,]$qstart,temp[f,]$qend))
    }
    mycounts<-table(target_seq2)-1 #subtract 1 because we've appended to the original sequence, which is necessary for calculating missingness
    paralog_percent<-length(mycounts[mycounts>1])/qlenuniq
    # weighted_paralog_percent<-(length(mycounts[mycounts>1])*mean(mycounts[mycounts>1])/qlenuniq)
    paralogy_index<-mean(mycounts)
    paralogy_index_ignoreMissing<-mean(mycounts[mycounts>0])    # this can be interpreted as copy number for recovered sequences
    full_percent<-length(mycounts[mycounts>0])/qlenuniq
    unique_percent<-length(mycounts[mycounts==1])/qlenuniq
    missing_percent<-length(mycounts[mycounts==0])/qlenuniq
    paralog_percent_ignoreMissing<-paralog_percent/full_percent
    coverage_summary_chromosome_UNaware<-rbind(coverage_summary_chromosome_UNaware,
                                               data.frame(qseqid=unique(temp$qseqid),
                                                          paralogy_index=round(paralogy_index,2),
                                                          paralogy_index_ignoreMissing=round(paralogy_index_ignoreMissing,2),
                                                          paralog_percent=round_perc(paralog_percent),
                                                          paralog_percent_ignoreMissing=round_perc(paralog_percent_ignoreMissing),
                                                          full_percent=round_perc(full_percent),
                                                          unique_percent=round_perc(unique_percent),
                                                          missing_percent=round_perc(missing_percent)))
  }
  if(doCovPerChrom) {
    # coverage_summary_chromosome_aware<-data.frame(cbind(coverage_summary_chromosome_aware[,1:4],apply(coverage_summary_chromosome_aware[,5:ncol(coverage_summary_chromosome_aware)],2,function(x) round(x,3)*100)))
    # coverage_summary_chromosome_aware<-data.frame(cbind(coverage_summary_chromosome_aware[,1:4],apply(coverage_summary_chromosome_aware[,grep("percent",names(coverage_summary_chromosome_aware))],2,function(x) round(x,3)*100)))
    # coverage_summary_chromosome_UNaware<-data.frame(cbind(coverage_summary_chromosome_UNaware[,1:3],apply(coverage_summary_chromosome_UNaware[,grep("percent",names(coverage_summary_chromosome_UNaware))],2,function(x) round(x,3)*100)))
    
    coverage_summary_nchroms<-coverage_summary_chromosome_aware %>%
      group_by(qseqid) %>%
      summarise(n_chroms=n(),
                .groups = "keep")
    coverage_summary_chromosome_UNaware<-left_join(coverage_summary_chromosome_UNaware,coverage_summary_nchroms,"qseqid")
    return(list(coverage_summary_chromosome_aware=coverage_summary_chromosome_aware,coverage_summary_chromosome_unaware=coverage_summary_chromosome_UNaware))
  }else{
    ## Note when doCovPerChrom=FALSE you don't get the n_chroms field
    # coverage_summary_chromosome_UNaware<-data.frame(cbind(coverage_summary_chromosome_UNaware[,1:3],apply(coverage_summary_chromosome_UNaware[,4:ncol(coverage_summary_chromosome_UNaware)],2,function(x) round(x,3)*100)))
    return(list(coverage_summary_chromosome_unaware=coverage_summary_chromosome_UNaware))
  }
}

FindIntrons<-function(data,max_intron_length=10000,max_intron_percent=100){
  introns_out<-data.frame()
  introns_flag_out<-data.frame()
  for(i in unique(data$qseqid)){
    
    temp<-data %>% 
      filter(qseqid==i)
    
    introns <-  temp %>% group_by(qseqid,sseqid)%>% 
      mutate(p_sstart=ifelse(sstart<send,sstart,send),
             p_send=ifelse(p_sstart==sstart,send,sstart),
             orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
      arrange(p_sstart,.by_group = TRUE) %>%
      summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
                intron_length=shift_sstart-shift_send,
                orientation=c(orientation,NA),
                intron_position_on_target=c(qend,NA),
                supercontig_length=max(p_send)-min(p_sstart),
                .groups = "keep")%>%
      na.exclude() %>%
      select(qseqid,sseqid,intron_length,intron_position_on_target,supercontig_length)
    ## write to data.frame
    introns_out<-rbind(introns_out,introns)
    
    introns_flag<-introns %>%
      summarise(longest_intron_length=max(intron_length),
                total_intron_length=sum(intron_length),
                largest_intron_percent=max(sum(intron_length)/unique(supercontig_length))*100,
                exceeds_max_intron_length=ifelse(any(intron_length>=max_intron_length),TRUE,FALSE),
                exceeds_max_intron_percent=ifelse(largest_intron_percent>=max_intron_percent,TRUE,FALSE),
                summed_supercontig_length=unique(supercontig_length),
                .groups="keep") %>%
      ungroup() %>%
      select(qseqid,sseqid,summed_supercontig_length,
             longest_intron_length,total_intron_length,exceeds_max_intron_length,
             largest_intron_percent,exceeds_max_intron_percent) %>%
      unique()
    
    introns_flag_out<-rbind(introns_flag_out,introns_flag)
  }
  return(list(intron_details=introns_out,targets_with_intron_flags=introns_flag_out))
}

PlotTargets<-function(data,output_prefix,min_display_intron,max_display_intron){
  pdf(paste0(output_prefix,"_plots.pdf"),width=15)
  cat("Making synteny plots. It's kinda hard work. Please bear with me ^_^ \n")
  pb <- progress_bar$new(total = length(unique(data$qseqid)),
                         clear=F,
                         format = "[:bar] :percent; elapsed :elapsed")
  for(i in unique(data$qseqid)){
    pb$tick()
    # cat("Plotting",i,"\n")
    temp<-data %>% 
      filter(qseqid==i)
    # don't plot if it won't fit in the plot area
    if(nrow(temp)<=72){
      introns <-  temp %>% group_by(qseqid,sseqid)%>% 
        mutate(p_sstart=ifelse(sstart<send,sstart,send),
               p_send=ifelse(p_sstart==sstart,send,sstart),
               orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
        arrange(p_sstart,.by_group = TRUE) %>%
        summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
                  intron=(shift_sstart-shift_send),
                  qpos=c(NA,qend),shift_qstart=c(qstart,NA),
                  orientation=c(orientation,NA),
                  .groups = "keep")%>%
        filter(intron>min_display_intron,intron<max_display_intron)%>%
        na.exclude()
      
      p<-ggplot(temp)+
        geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=pident),
                     show.legend = T,size=2.5)+
        scale_colour_viridis_c("Sequence\nsimilarity")+
        facet_wrap(~sseqid,scales = "free")+
        xlim(c(0,unique(temp$qlen)))+
        theme_bw()+
        labs(x="Target position",y="Genome position")+
        ggtitle(i)+
        coord_flip()
      if(nrow(introns)>0){
        p<-p+geom_rect(data=introns,
                       aes(xmin=-Inf,xmax=Inf,ymin=shift_sstart,ymax=shift_send),
                       colour="darkgreen",alpha=0.05,size=NA)+
          geom_label_repel(data=introns,
                           aes(x=qpos,
                               y=shift_sstart-(intron/2),
                               label=intron),fill="white",
                           force_pull = 0.1,
                           size=2,
                           min.segment.length=unit(5,"mm"),
                           box.padding = 0.05,label.padding=0.07,
                           label.r=0.01,seed=1234,
                           max.overlaps=15)
      }
      suppressWarnings(print(p))
    }else{
      p<-ggplot()+ggtitle(paste0(i,": Too many matches to fit in plot area!"))
      suppressWarnings(print(p))
    }
  }
  dev.off()
}

enumerate<-function(x){seq(min(x),max(x),by=1)}
full_enumeration<-function(x,a,b){do.call(c,sapply(1:nrow(x),function(k){enumerate(x[k,c(a,b)])},simplify = F))}
full_containment<-function(x,a,b){mean(table(full_enumeration(x,a,b)))}

ThinBlastResult<-function(data){
  cat("Thinning blast result...\n")
  dat<-data
  dat$pair<-paste0(dat$qseqid,"_",dat$sseqid)
  new.dat<-data.frame()
  pb <- progress_bar$new(total = length(unique(dat$pair)),
                         clear = F,
                         show_after = 0,
                         format = "[:bar] :percent; elapsed :elapsed")
  # kombis_out<-list()
  for(i in unique(dat$pair)){
    pb$tick()
    dat.thin<-dat[dat$pair==i,] %>% arrange(desc(length))
    dat.thin.containment<-full_containment(dat.thin,"qstart","qend")
    
    if(nrow(dat.thin)>1 & dat.thin.containment>=1.1){
      kombis<-data.frame(t(combn(1:nrow(dat.thin),m=2)))
      kombis$containment<-sapply(1:nrow(kombis),function(g) full_containment(dat.thin[c(kombis[g,1],kombis[g,2]),],"qstart","qend"))
      kombis$pair<-i
      if(any(kombis$containment>=1.1)){
        kombis.todo<-kombis[kombis$containment>=1.1,]
        kombis.todo$worstscoring_index<-NA
        for(g in 1:nrow(kombis.todo)){
          my_pair<-c(kombis.todo[g,1],kombis.todo[g,2])
          worstscoring_index<-ifelse(length(unique(dat.thin[my_pair,"bitscore"]))==2,
                                                     yes=which(dat.thin[my_pair,"bitscore"]==min(dat.thin[my_pair,"bitscore"])),
                                                     no=ifelse(length(unique(dat.thin[my_pair,"pident"]))==2,
                                                               which(dat.thin[my_pair,"pident"]==min(dat.thin[my_pair,"pident"])),
                                                               which(dat.thin[my_pair,"length"]==min(dat.thin[my_pair,"length"]))
                                                               )
          )## if bitscores are different, remove lowest, otherwise (very unlikely) remove lowest pident, and if THOSE are identical, remove shortest, and if THOSE are identical, remove second
          worstscoring_index<-ifelse(length(worstscoring_index)!=1,worstscoring_index[2],worstscoring_index)
          kombis.todo[g,]$worstscoring_index<-worstscoring_index
        }
        kombis.todo$worstscoring<-sapply(1:nrow(kombis.todo),function(x) {kombis.todo[x,kombis.todo[x,"worstscoring_index"]]},simplify = "array")
        dat.thin<-dat.thin[-unique(as.numeric(kombis.todo$worstscoring)),]
        # dat.thin.containment<-full_containment(dat.thin,"qstart","qend")
      }
      # kombis_out[[i]]<-kombis.todo
    }
    new.dat<-rbind(new.dat,dat.thin)
  }
  cat("Original blast result had",nrow(dat),"rows. After thinning",
        nrow(new.dat),"rows (",round(nrow(new.dat)/nrow(dat)*100,1) ,"%) remain.\n")
  # return(list(thinned_blast_result=new.dat,
  #             combinations_list=kombis_out))
  return(new.dat)
}

CheckTargets<-function(blast_file,
                       min_pident,
                       min_fragment_length,
                       max_intron_length,
                       max_intron_percent,
                       output_prefix,
                       min_display_intron,
                       max_display_intron,
                       doPlots,
                       doIntronStats,
                       doCovPerChrom,
                       multicopyTarget,
                       genelist,
                       blast_type,
                       doThin){
                         
  # suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
  suppressMessages(suppressWarnings(require(dplyr,quietly=TRUE,warn.conflicts=FALSE)))
  suppressMessages(suppressWarnings(require(tidyr,quietly=TRUE,warn.conflicts=FALSE)))
  suppressMessages(suppressWarnings(require(progress,quietly=TRUE,warn.conflicts=FALSE)))

  cat("\nBeginning VetTargets_genome.R...\n")
  cat("Reading blast file:", blast_file,".\n")
  dat <- as_tibble(read.table(blast_file,header=T))
  if(doThin){
    if(!file.exists(paste0(output_prefix,"_",blast_type,"_Thinned_minPident",min_pident,"_minLength",min_fragment_length,".txt"))){
      cat("Now thinning BLAST result, which has",nrow(dat),"rows.\n")
      cat("BLAST type specified as",blast_type,".\n")
      if(blast_type=="tblastx"){
        cat("Multiplying length by 3 to get length in nucleotides rather than amino acids.\n")
        cat("Will remove hits <",min_fragment_length,"nucleotide base pairs long.\n")
        dat$length.nuc<-dat$length*3
      } else {
        dat$length.nuc<-dat$length
      }
      dat <- dat %>%
        filter(pident >= min_pident,
              length.nuc >= min_fragment_length)
      cat("After removing matches with pident <",min_pident,"and length <",min_fragment_length, "BLAST result has",nrow(dat),"rows.\n")
      cat("Now removing redundant BLAST hits.\n")
      dat <- ThinBlastResult(dat)
      
      write.table(dat,file = paste0(output_prefix,"_",blast_type,"_Thinned_minPident",min_pident,"_minLength",min_fragment_length,".txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
      cat("After removing redundant hits, BLAST result has",nrow(dat),"rows.\n")
      cat("This thinned blast result was written to",paste0(output_prefix,"_",blast_type,"_Thinned_minPident",min_pident,"_minLength",min_fragment_length,".txt"),".\n")
    } else {
      cat("Thinned BLAST result already exists. Will read it in instead of repeating the thinning procedure.\n")
      cat("Reading in",paste0(output_prefix,"_",blast_type,"_Thinned_minPident",min_pident,"_minLength",min_fragment_length,".txt"),".\n")
      dat<-as_tibble(read.table(paste0(output_prefix,"_",blast_type,"_Thinned_minPident",min_pident,"_minLength",min_fragment_length,".txt"),header = T))
    }
  }
  

  if(!is.null(genelist)){
    gl <- suppressMessages(scan(genelist,what="character"))
  }
  ## for cases where we have multiple copies of each gene in the target file, 
  # split the blast results by copy source and run analyses separately on each.
  if(multicopyTarget){
    dat <- dat %>% separate(col = qseqid, into = c("Source","qseqid"), sep = "-")
    sp <- dat$Source
    dat.list <- split(dat, sp)
    output_prefix_OG <- output_prefix
    for (i in names(dat.list)){
      cat("\nTargets are represented by multiple sources. Now working on targets from", i, "...\n")
      dat <- dat.list[[i]]
      if(!is.null(genelist)){
        dat <- dat[dat$qseqid %in% gl,]
        if(nrow(dat)==0) stop("After filtering out genes not present in the provided gene list, nothing remains! Check the names match.\n")
      }
      ## put results for each copy in separate directory
      if(!dir.exists(i)) {
        dir.create(i)
      }
      output_prefix <- paste0(i,"/",output_prefix_OG)
      if(file.exists(paste0(output_prefix,"_CoverageStats_AcrossChromosomes.txt"))){
        cat(paste0(output_prefix,"_CoverageStats_AcrossChromosomes.txt"),"exists. Skipping.\n")
      }else{
        cov_stats<-GetCoverageStats(dat,doCovPerChrom)
        if(doCovPerChrom){ 
          write.table(cov_stats$coverage_summary_chromosome_aware,paste0(output_prefix,"_CoverageStats_PerChromosome.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" ) 
        }
        write.table(cov_stats$coverage_summary_chromosome_unaware,paste0(output_prefix,"_CoverageStats_AcrossChromosomes.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
        
        if(doIntronStats){
          intron_stats<-suppressWarnings(FindIntrons(data=dat,max_intron_length,max_intron_percent))
          write.table(intron_stats$intron_details,paste0(output_prefix,"_IntronStats.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
          write.table(intron_stats$targets_with_intron_flags,paste0(output_prefix,"_IntronFlags.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
        }
        
        if(doPlots){
          suppressMessages(suppressWarnings(require(ggrepel,quietly=TRUE,warn.conflicts=FALSE)))
          PlotTargets(dat,output_prefix,min_display_intron,max_display_intron)
        }
      }
    }
    
  } else {
    if(!is.null(genelist)){
      dat <- dat[dat$qseqid %in% gl,]
      if(nrow(dat)==0) stop("After filtering out genes not in the provided gene file, nothing remains! Check that the names match.\n")
    }
    cov_stats<-GetCoverageStats(dat,doCovPerChrom)
    if(doCovPerChrom){ write.table(cov_stats$coverage_summary_chromosome_aware,paste0(output_prefix,"_CoverageStats_PerChromosome.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" ) }
    write.table(cov_stats$coverage_summary_chromosome_unaware,paste0(output_prefix,"_CoverageStats_AcrossChromosomes.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
    
    if(doIntronStats){
      intron_stats<-suppressWarnings(FindIntrons(data=dat,max_intron_length,max_intron_percent))
      write.table(intron_stats$intron_details,paste0(output_prefix,"_IntronStats.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
      write.table(intron_stats$targets_with_intron_flags,paste0(output_prefix,"_IntronFlags.txt"),quote = F,row.names = F,col.names = TRUE, sep = "\t" )
    }
    
    if(doPlots){
      suppressMessages(suppressWarnings(require(ggrepel,quietly=TRUE,warn.conflicts=FALSE)))
      PlotTargets(dat,output_prefix,min_display_intron,max_display_intron)
    }
  }
}

#### DONE DEFINING FUNCTIONS
suppressMessages(suppressWarnings(require(optparse)))

p <- OptionParser(usage="This script will take a tabular blast result (-outfmt 6) **WITH A HEADER LINE** with your target fasta as the query and (draft) genome as the subject\n
  and do the following:\n
   1. Find potential paralogs by identifying genes with good matches from more than one contig/chromosome segments in the reference genome\n
   2. Find potentially missing genes by identifying genes with no or few good matches in the reference genome\n
   3. if --doIntronStats=TRUE, find genes spanning huge introns\n
  then\n
   4. if --doPlots=TRUE, create sets of simple segment plots to allow for visual inspection of blast alignments of genes failing and passing the above checks\n
   5. Output files listing names of targets passing or failing each check, and a file of 'clean' targets passing all checks\n
  Run using Rscript, e.g.\n
  Rscript VetTargets_genome.R --blast_file blastn_targets_to_genome.txt --min_fragment_length 50 --min_pident 80 --max_intron_length 10000 --max_intron_percent 60 --output_prefix test\n")
# Add a positional argument
# Required
p <- add_option(p, c("-b","--blast_file"), help="<Required: tab-delimited blast result, target=query, genome=subject>",type="character")
p <- add_option(p, c("-f","--min_fragment_length"), help="<Required: minimum length of a blast hit, ignore anything shorter>",type="numeric")
p <- add_option(p, c("-p","--min_pident"), help="<Required: minimum % identity of a blast hit, ignore anything lower>",type="numeric")
p <- add_option(p, c("-i","--max_intron_length"), help="<Required: length threshold of the largest intron in a gene; flag any gene exceeding>",type="numeric")
p <- add_option(p, c("-I","--max_intron_percent"), help="<Required: percentage of a supercontig's length consisting of introns; flag any gene exceeding>",type="numeric")
# Optional
p <- add_option(p, c("-o","--output_prefix"), help="<prefix to name results files; defaults to VetTargets_genome_results>",type="character",default = "VetTargets_genome_results")
p <- add_option(p, c("-d","--min_display_intron"), help="<intron length above which to annotate on plots; default 1kb>",type="numeric",default=1000)
p <- add_option(p, c("-D","--max_display_intron"), help="<don't annotate introns longer than this; default 1Mb>",type="numeric",default=1e6)
p <- add_option(p, c("-P","--doPlots"), help="<Make plots? Default=TRUE>",type="logical",default=TRUE)
p <- add_option(p, c("-S","--doIntronStats"), help="<Calculate intron stats? Default=TRUE>",type="logical",default=TRUE)
p <- add_option(p, c("-C","--doCovPerChrom"), help="<Calculate per-chromosome coverage stats (in addition to across-chromosome)? Default=TRUE>",type="logical",default=TRUE)
p <- add_option(p, c("-M","--multicopyTarget"), help="<does target file contain multiple copies per gene (TRUE or FALSE)? If TRUE, gene names must follow HybPiper convention, E.g. Artocarpus-gene001 and Morus-gene001 are the same gene. Default=FALSE>",type="logical",default=FALSE)
p <- add_option(p, c("-g","--genelist"), help="<file listing genes to process, excluding any others>",type="character",default=NULL)
p <- add_option(p, c("-B","--blast_type"), help="<blastn or tblastx? Default=blastn>",type="character",default="blastn")
p <- add_option(p, c("-T","--doThin"), help="<whether to thin blast results>",type="logical",default=TRUE)

# parse
args<-parse_args(p)
# print(args)
if(is.null(args$blast_file)){ 
  print_help(p)
  }else{
  ## RUN
  try(CheckTargets(blast_file = args$blast_file,
                   min_pident = args$min_pident,
                   min_fragment_length = args$min_fragment_length,
                   max_intron_length = args$max_intron_length,
                   max_intron_percent = args$max_intron_percent,
                   output_prefix = args$output_prefix,
                   min_display_intron = args$min_display_intron,
                   max_display_intron = args$max_display_intron,
                   doPlots = args$doPlots,
                   doIntronStats = args$doIntronStats,
                   doCovPerChrom = args$doCovPerChrom,
                   multicopyTarget = args$multicopyTarget,
                   genelist = args$genelist,
                   blast_type = args$blast_type,
                   doThin = args$doThin))
  cat("VetTargets_genome.R is done!\n\n#############################################")
}

