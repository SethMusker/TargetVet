
collate<-function(samples,directory){
    
    samples<-scan(samples,what="character")
    
    out<-data.frame()
    for(i in samples){
        temp<-read.table(paste0(directory,"/",i,"_CoverageStats_AcrossChromosomes.txt"),header=T)
        temp$Sample<-i
        out<-rbind(out,temp)
    }
    out<-as_tibble(out)
    out$Sample<-as.factor(out$Sample)
    out$qseqid<-as.factor(out$qseqid)
    
    out_meanSort<-out %>% mutate(qseqid=fct_reorder(qseqid,paralog_percent_ignoreMissing,mean)) %>%
        arrange(paralog_percent_ignoreMissing)
    
    #get the order and do segmented
    out_meanSort$orderMean<-sapply(out_meanSort$qseqid,function(x) {
        y<-which(x==levels(out_meanSort$qseqid))
        return(y)})
    lm.mean<-lm(paralog_percent_ignoreMissing~orderMean+0,data=out_meanSort) # force intercept=0
    seg.mean<-segmented(lm.mean,npsi=1)
    out_meanSort$seg.pred<-predict(seg.mean,data.frame(orderMean=out_meanSort$orderMean))
    
    #####################
    ### BEGIN PLOTTING ##
    #####################
    out_meanSort %>%
        ggplot(aes(x=qseqid,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_point(aes(colour=missing_percent),alpha=0.5)+
        scale_colour_gradientn(colours=c("blue","grey","red"))+
        theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5))+
        geom_line(aes(x=orderMean,y=seg.pred))+
        geom_vline(xintercept = seg.mean$psi[2],lty=2)+
        labs(y="Paralog Rate (%)")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_segmented.pdf"),width=40,height=20,units = "cm")
    
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        ylim(c(-1,1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("All target loci")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_segmented_RESIDUALS_boxplot.pdf"),width=30,units = "cm")
    
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        ylim(c(-1,1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_segmented_RESIDUALS_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")
    #####################
    ### DONE PLOTTING ###
    #####################
    
    out_meanSort_summary<-out_meanSort %>% group_by(qseqid) %>%
        summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                    mean_missing_percent=round(mean(missing_percent),1),
                    mean_paralogy_index=round(mean(paralogy_index),1),
                    diagnosis=ifelse(unique(orderMean)<=seg.mean$psi[2],"Paralog","Single"),
                    .groups="keep") 
    
    write.table(out_meanSort_summary,paste0(directory,"/paralogy_summary.txt"),
                quote=F,row.names=F)
    
    out_meanSort_summary %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
        write.table(paste0(directory,"/paralog_list.txt"),
                    quote=F,row.names=F,col.names=F)
    out_meanSort_summary %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
        write.table(paste0(directory,"/singleCopy_list.txt"),
                    quote=F,row.names=F,col.names=F)
}

suppressMessages(suppressWarnings(require(optparse,quietly=TRUE,warn.conflicts=FALSE)))

p <- OptionParser(usage=" This script will tale the CoverageStats output of VetTargets_genome.R for \n
                        many samples and detect paralogs in a targeted set of genes ")
# Add a positional argument
p <- add_option(p, c("-s","--samples"), help="<Required: list of sample names>",type="character")
p <- add_option(p, c("-d","--directory"), help="<Required: directory with output from VetTargets>",type="character")
# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(segmented,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))

collate(samples=args$samples,
        directory=args$directory)

















