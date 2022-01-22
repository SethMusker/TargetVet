#!/bin/Rscript
VetTargets_WGS<-function(cov,target,out,method,SD_cutoff,max_median=NULL,min_median=NULL){
  if(method == "sd"){method<-"SD"}
  cov <- unlist(strsplit(cov, ","))
  incov<-list()
  incov_summary<-list()
  chosen<-list()
  highcov<-list()
  lowcov<-list()
  for(i in cov){
      incov[[i]]<-read.table(i,stringsAsFactors = T,header=F)
      names(incov[[i]])<-c("ID","pos","coverage")
        incov_summary[[i]]<-incov[[i]] %>%
            group_by(ID) %>%
            summarise(median_cov=median(coverage),
                        sd_cov=round(sd(coverage),1)) %>%
            arrange(dplyr::desc(median_cov)) %>%
            mutate(ID = fct_reorder(ID, dplyr::desc(median_cov)))
        incov_summary[[i]]$ord<-c(1:nrow(incov_summary[[i]]))
        write.table(incov_summary[[i]],paste0(i,"_summary.txt"),row.names=F,quote=F)
        if(method=="breakpoint"){
          ## breakpoint analysis
          myseg<-segmented(lm(median_cov~ord,data=incov_summary[[i]]),npsi=2)
          breakpoint_high_mc<-as.numeric(incov_summary[[i]][floor(myseg$psi[1,2]),"median_cov"])
          breakpoint_low_mc<-as.numeric(incov_summary[[i]][ceiling(myseg$psi[2,2]),"median_cov"])
          incov_summary[[i]]$segment_line<-predict(myseg)
          incov_summary[[i]]$segment_line[incov_summary[[i]]$segment_line < 0]<-0
        }
        if(method=="SD"){
          breakpoint_high_mc<-  floor(mean(incov_summary[[i]]$median_cov) + SD_cutoff*sd(incov_summary[[i]]$median_cov))
          breakpoint_low_mc<- ceiling(mean(incov_summary[[i]]$median_cov) - SD_cutoff*sd(incov_summary[[i]]$median_cov))
        }
        if(method=="manual"){
          breakpoint_high_mc<-max_median
          breakpoint_low_mc<-min_median
        }

        cat("Cut-off values for", i,": \nupper median coverage =",breakpoint_high_mc,"\nlower median coverage =",breakpoint_low_mc,"\n")
            
        chosen[[i]]<-  incov_summary[[i]] %>% dplyr::filter(median_cov >= breakpoint_low_mc, median_cov <= breakpoint_high_mc) %>% select(ID) %>% as.matrix() %>% as.vector()
        highcov[[i]]<-incov_summary[[i]] %>% dplyr::filter(median_cov > breakpoint_high_mc) %>% select(ID) %>% as.matrix() %>% as.vector()
        lowcov[[i]]<- incov_summary[[i]] %>% dplyr::filter(median_cov < breakpoint_low_mc) %>% select(ID) %>% as.matrix() %>% as.vector()
        
        ## make plots

        myplot<-ggplot(incov_summary[[i]])+
          geom_col(aes(x=ID,y=median_cov,fill=sd_cov))+
          scale_fill_viridis_c("Coverage SD")+
          theme_bw()+
          theme(axis.text.x = element_text(angle=90,size=rel(0.5)),
                panel.grid = element_blank())+
          scale_y_continuous(n.breaks = 25)+
          geom_hline(yintercept = breakpoint_high_mc,colour="red")+
          geom_hline(yintercept = breakpoint_low_mc,colour="red")+
          ylim(c(0,max(incov_summary[[i]]$median_cov)))+
          labs(y="Median Coverage")+
          ggtitle(paste0("Method: ",method))
        
        if(method=="breakpoint"){
          myplot<-myplot+
            geom_line(aes(x=ord,y=segment_line),size=1,colour="skyblue",alpha=0.8)
            
        }
        ggsave(paste0("MedianCoverage_",out,"_method_",method,".pdf"),plot=myplot,width=30,height=20,units="cm")

    }
    
  if(length(cov)>1){
    chosen_intersect<-Reduce(intersect, chosen)
    chosen_outersect<-incov_summary[[i]]$ID[! incov_summary[[i]]$ID %in% chosen_intersect]
    highcov_intersect<-Reduce(intersect, highcov)
    lowcov_intersect<-Reduce(intersect, lowcov)
  }else{
    chosen_intersect<-chosen[[1]]
    chosen_outersect<-incov_summary[[1]]$ID[! incov_summary[[1]]$ID %in% chosen_intersect]
    highcov_intersect<- highcov[[1]]
    lowcov_intersect<- lowcov[[1]]
  }
    cat("Filtering kept",length(chosen_intersect),"out of",length(incov_summary[[i]]$ID),"targets")
    write.table(data.frame(ID=chosen_intersect),paste0(out,".kept"),col.names=F,row.names=F,quote=F)
    write.table(data.frame(ID=chosen_outersect),paste0(out,".removed"),col.names=F,row.names=F,quote=F)
    write.table(data.frame(ID=highcov_intersect),paste0(out,".highcoverage"),col.names=F,row.names=F,quote=F)
    write.table(data.frame(ID=lowcov_intersect),paste0(out,".lowcoverage"),col.names=F,row.names=F,quote=F)
    
    mytarget<-readDNAStringSet(target)
    myout_target<-mytarget[chosen_intersect]
    writeXStringSet(myout_target,paste0(out,"_middle_coverage.fasta"))
    
    myout_target<-mytarget[highcov_intersect]
    writeXStringSet(myout_target,paste0(out,"_high_coverage.fasta"))
    
    myout_target<-mytarget[lowcov_intersect]
    writeXStringSet(myout_target,paste0(out,"_low_coverage.fasta"))
}

suppressMessages(suppressWarnings(require(optparse)))

p <- OptionParser(usage="Vet your targets by identifying loci with low/null, roughly expected, and greater than expected median coverage.\n
                Requires the following packages: tidyverse, Biostrings, segmented, optparse\n
                Run using Rscript, e.g.\n
                Rscript VetTargets.R --cov 1st_BAM_to_myTARGETS.coverage,2nd_BAM_to_myTARGETS.coverage --target myTARGETS.fasta --out myTARGETS_filtered")
# Add a positional argument
p <- add_option(p, c("-c","--cov"), help="REQUIRED coverage file(s) output by map_WGS_to_targets.sh or made by you\n(if so, make sure it has all sites included - e.g. by using samtools depth -a)\nif you have coverage files for >1 species/sample, just write them all in here, comma separated (NO SPACES)")
p <- add_option(p, c("-t","--target"), help="REQUIRED target file used as the reference for the coverage calculation")
p <- add_option(p, c("-o","--out"), help="REQUIRED: Prefix of filtered fasta files to write")
p <- add_option(p, c("-m","--method"), help="REQUIRED: method to determine high and low median coverage values used to classify genes:\n one of 'breakpoint' (default), 'SD', 'manual'",default="breakpoint")
p <- add_option(p, c("-s","--SD_cutoff"), help="if method=SD or sd, 'good' genes will fall within mean(median coverage) +- SD_cutoff*sd(median coverage). Default = 1",default=1)
p <- add_option(p, c("-u","--max_median"), help="REQUIRED IF method=manual, 'good' genes will fall between max_median and min_median")
p <- add_option(p, c("-l","--min_median"), help="REQUIRED IF method=manual. See above.")

# parse
args<-parse_args(p)

# suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dplyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(forcats,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(Rsamtools,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(Biostrings,quietly =TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(segmented,quietly =TRUE,warn.conflicts=FALSE)))

VetTargets_WGS(cov=args$cov,
  target=args$target,
  out=args$out,
  method=args$method,
  SD_cutoff=args$SD_cutoff,
  max_median=args$max_median,
  min_median=args$min_median)




