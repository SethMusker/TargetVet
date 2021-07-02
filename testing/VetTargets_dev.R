#!/bin/Rscript
suppressMessages(suppressWarnings(require(optparse)))

filterTargets<-function(cov,target,filterParalogs,out){
  cov <- unlist(strsplit(cov, ","))
  incov<-list()
  incov_summary<-list()
  mychosen<-list()
  for(i in cov){
      incov[[i]]<-read.table(i,stringsAsFactors = T,header=F)
      names(incov[[i]])<-c("ID","pos","coverage")
      ## breakpoint analysis
        incov_summary[[i]]<-incov[[i]] %>%
            group_by(ID) %>%
            summarise(median_cov=median(coverage),
                        sd_cov=round(sd(coverage),1)) %>%
            arrange(dplyr::desc(median_cov)) %>%
            mutate(ID = fct_reorder(ID, dplyr::desc(median_cov)))
        incov_summary[[i]]$ord<-c(1:nrow(incov_summary[[i]]))
        myseg<-segmented(lm(median_cov~ord,data=incov_summary[[i]]),npsi=2)
        breakpoint_high_mc<-as.numeric(incov_summary[[i]][floor(myseg$psi[1,2]),"median_cov"])
        breakpoint_low_mc<-as.numeric(incov_summary[[i]][ceiling(myseg$psi[2,2]),"median_cov"])
        cat(i,"\nupper=",breakpoint_high_mc,"\nlower=",breakpoint_low_mc,"\n")
        if(filterParalogs==0){
            mychosen[[i]]<-incov_summary[[i]] %>% dplyr::filter(median_cov >= breakpoint_low_mc) %>% select(ID) %>% as.matrix() %>% as.vector()
        } else {
            mychosen[[i]]<-incov_summary[[i]] %>% dplyr::filter(median_cov >= breakpoint_low_mc, median_cov <= breakpoint_high_mc) %>% select(ID) %>% as.matrix() %>% as.vector()
        }
        ## make plots
        incov_summary[[i]]$segment_line<-predict(myseg)
        incov_summary[[i]]$segment_line[incov_summary[[i]]$segment_line < 0]<-0
        myplot<-ggplot(incov_summary[[i]])+
            geom_col(aes(x=ID,y=median_cov,fill=sd_cov))+
            scale_fill_viridis_c("Coverage SD")+
            theme(axis.text.x  = element_text(angle=90,size=5.5))+
            scale_y_continuous(n.breaks = 25)+
            geom_hline(yintercept = breakpoint_high_mc,colour="red")+
            geom_hline(yintercept = breakpoint_low_mc,colour="red")+
            geom_line(aes(x=ord,y=segment_line),size=1,colour="green",alpha=0.8)+
            ylim(c(0,max(incov_summary[[i]]$median_cov)))+
            ylab("Median Coverage")
        ggsave(paste0(i,".pdf"),plot=myplot,width=25,height=20,units="cm")

    }
    mychosen_intersect<-Reduce(intersect, mychosen)
    mychosen_outersect<-incov_summary[[i]]$ID[! incov_summary[[i]]$ID %in% mychosen_intersect]
    cat("Filtering kept",length(mychosen_intersect),"out of",length(incov_summary[[i]]$ID),"targets")
    write.table(data.frame(ID=mychosen_intersect),paste0(out,".kept"),col.names=F,row.names=F,quote=F)
    write.table(data.frame(ID=mychosen_outersect),paste0(out,".removed"),col.names=F,row.names=F,quote=F)
    myout_target<-readDNAStringSet(target)[mychosen_intersect]
    writeXStringSet(myout_target,out)
}


p <- OptionParser(usage="Vet your targets by removing loci with low/null median coverage and (optionally) loci with high coverage (potential paralogs).\n
                Requires the following packages: tidyverse, Biostrings, segmented, optparse\n
                Run using Rscript, e.g.\n
                Rscript VetTargets.R --cov 1st_BAM_to_myTARGETS.coverage,2nd_BAM_to_myTARGETS.coverage --target myTARGETS.fasta --filterParalogs 1 --out myTARGETS_filtered.fasta")
# Add a positional argument
p <- add_option(p, c("-c","--cov"), help="coverage file(s) output by map_WGS_to_targets.sh or made by you\n(if so, make sure it has all sites included - e.g. by using samtools depth -a)\nif you have coverage files for >1 species/sample, just write them all in here, comma separated (NO SPACES)")
p <- add_option(p, c("-t","--target"), help="your target file used as the reference for the coverage calculation")
p <- add_option(p, c("-f","--filterParalogs"), help="filter putative paralogs? 0:no, 1:yes",default=1,type="numeric")
p <- add_option(p, c("-o","--out"), help="name of filtered fasta file (including extension) to write")

# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(tidyverse,quietly =TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(Biostrings,quietly =TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(segmented,quietly =TRUE,warn.conflicts=FALSE)))
filterTargets(cov=args$cov,target=args$target,filterParalogs=args$filterParalogs,out=args$out)




