
collate<-function(samples,directory){
    
    samples<-scan(samples,what="character")
    
    out<-data.frame()
    for(i in samples){
        temp<-read.table(paste0(directory,"/",i,"_CoverageStats_AcrossChromosomes.txt"),header=T)
        # if("n_chroms" %in% names(temp)){ temp<-temp[,-which(names(temp)=="n_chroms")]}
        temp$Sample<-i
        out<-dplyr::bind_rows(out,temp)
    }
    out$Sample<-as.factor(out$Sample)
    out$qseqid<-as.factor(out$qseqid)
    out<-as_tibble(out)
    
    out_meanSort<-out %>%
        mutate(qseqid=fct_reorder(qseqid,paralog_percent_ignoreMissing,mean)) %>%
        arrange(paralog_percent_ignoreMissing)
    
    #get the order and do segmented
    out_meanSort$orderMean<-sapply(out_meanSort$qseqid,function(x) {
        y<-which(levels(out_meanSort$qseqid)==x)
        return(y)})
    lm.mean<-lm(paralog_percent_ignoreMissing~orderMean+0,data=out_meanSort) # force intercept=0
    seg.mean<-segmented(lm.mean,npsi=1)
    out_meanSort$seg.pred<-predict(seg.mean,newdata=data.frame(orderMean=out_meanSort$orderMean))
    out_meanSort$resid<-out_meanSort$paralog_percent_ignoreMissing - out_meanSort$seg.pred
    out_meanSort$resid_sq<-out_meanSort$resid^2
    
    #####################
    ### BEGIN PLOTTING ##
    #####################
    out_meanSort %>%
        ggplot(aes(x=qseqid,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_point(aes(colour=missing_percent),alpha=0.5)+
        scale_colour_gradientn("Missingness (%)",colours=c("blue","grey","red"))+
        theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5))+
        geom_line(aes(x=orderMean,y=seg.pred))+
        geom_vline(xintercept = seg.mean$psi[2],lty=2)+
        labs(y="Paralog Rate (%)",x="Target")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_breakpoint.pdf"),width=40,height=20,units = "cm")

## Per-sample plots
   #plot missingness
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(directory,"/missing_percent_perSample.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(directory,"/missing_percent_perSample_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot paralogy
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_perSample.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_perSample_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot residuals
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("All target loci")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_breakpoint_RESIDUALS_boxplot.pdf"),width=30,units = "cm")
    
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_breakpoint_RESIDUALS_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")
    #####################
    ### DONE PLOTTING ###
    #####################
    
    out_meanSort_summary<-out_meanSort %>% group_by(qseqid) %>%
        summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                    mean_missing_percent=round(mean(missing_percent),1),
                    mean_paralogy_index=round(mean(paralogy_index),1),
                    diagnosis=ifelse(unique(orderMean)<=seg.mean$psi[2],"Single","Paralog"),
                    .groups="keep")
    
    write.table(out_meanSort_summary,paste0(directory,"/paralogy_summary.txt"),
                quote=F,row.names=F)
    out_meanSort %>% group_by(qseqid) %>% arrange(paralog_percent_ignoreMissing,.by_group=TRUE) %>%
    write.table(paste0(directory,"/paralogy_details.txt"),
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

















