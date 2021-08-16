
collate<-function(samples,directory,outdir,force){
    
    if(dir.exists(outdir)){
        if(!force){
            stop("Specified output directory exists.\nTerminating. Use -f TRUE to overwrite.\n")
        }else{
            cat("Specified output directory exists.\n")
            cat("Will overwrite. \n(If it fails it's because the output file is open on your computer somewhere. Close and rerun.)\n")
        }
    }else{
        cat("Specified output directory exists not. Will attempt to create it.\n")
        dir.create(outdir)
    }
    
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
    lm.mean<-lm(paralog_percent_ignoreMissing~orderMean,data=out_meanSort) # DON'T force intercept=0
    seg.mean<-segmented(lm.mean,npsi=1)
    out_meanSort$seg.pred<-predict(seg.mean,newdata=data.frame(orderMean=out_meanSort$orderMean))
    out_meanSort$resid<-out_meanSort$paralog_percent_ignoreMissing - out_meanSort$seg.pred
    out_meanSort$resid_sq<-out_meanSort$resid^2
    
    #####################
    ### BEGIN PLOTTING ##
    #####################
        
    ##plot breakpoint analysis
    out_meanSort %>%
        ggplot(aes(x=qseqid,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_point(aes(colour=missing_percent),alpha=0.5)+
        scale_colour_viridis_c("Missingness (%)",option="D")+
        theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5))+
        geom_line(aes(x=orderMean,y=seg.pred))+
        geom_vline(xintercept = seg.mean$psi[2],lty=2)+
        labs(y="Paralog Rate (%)",x="Target")
    ggsave(paste0(directory,"/paralog_percent_ignoreMissing_meanSort_breakpoint.pdf"),width=40,height=20,units = "cm")
    
    ## HEATMAPS
    # paralogy
    ggplot(out_meanSort)+
        geom_raster(aes(x=qseqid,y=Sample,fill=paralog_percent_ignoreMissing))+
        scale_fill_viridis_c("Paralogy (%)",option="D")+
        labs(x="Target")+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3.5))
    ggsave(paste0(outdir,"/","paralogy_heatmap.pdf"),width=28,height=20,units="cm")
    # clustered heatmap (using heatmaply)
    out_mat_samp<-out_meanSort %>% select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = paralog_percent_ignoreMissing, values_fill = NA)
    out_mat_mat<-as.data.frame(out_mat_samp)
    rm(out_mat_samp)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    ggh<-ggheatmap(as.matrix(out_mat_mat),
                   k_row = 2, k_col = 2, 
                   seriate = "mean", 
                   fontsize_row = 6, fontsize_col = 5, 
                   column_text_angle = 90,  
                   na.value = "white",
                   key.title="Paralogy (%)")
    ggsave(paste0(outdir,"/","paralogy_heatmap_clustered.pdf"),width=28,height=20,units="cm",plot=ggh)
    
    #missingness
    out_meanSort %>% mutate(qseqid=fct_reorder(qseqid,missing_percent,mean)) %>% 
        ggplot()+
        geom_raster(aes(x=qseqid,y=Sample,fill=missing_percent))+
        scale_fill_viridis_c("Missingness (%)",option="D")+
        labs(x="Target")+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3.5))
    ggsave(paste0(outdir,"/","missingness_heatmap.pdf"),width=28,height=20,units="cm")
    # clustered heatmap (using heatmaply)
    out_mat_samp<-out_meanSort %>% select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = missing_percent, values_fill = NA)
    out_mat_mat<-as.data.frame(out_mat_samp)
    rm(out_mat_samp)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    ggh<-ggheatmap(as.matrix(out_mat_mat),
                   k_row = 2, k_col = 2, 
                   seriate = "mean", 
                   fontsize_row = 6, fontsize_col = 5, 
                   column_text_angle = 90,  
                   na.value = "white",
                   key.title="Missingness (%)")
    ggsave(paste0(outdir,"/","missingness_heatmap_clustered.pdf"),width=28,height=20,units="cm",plot=ggh)
    
    
## Per-sample plots
   #plot missingness
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/missing_percent_perSample.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/missing_percent_perSample_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot paralogy
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_perSample.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_perSample_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot residuals
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_meanSort_breakpoint_RESIDUALS_boxplot.pdf"),width=30,units = "cm")
    
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_meanSort_breakpoint_RESIDUALS_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")
    #####################
    ### DONE PLOTTING ###
    #####################
    
    out_meanSort_summary<-out_meanSort %>% group_by(qseqid) %>%
        summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                    mean_missing_percent=round(mean(missing_percent),1),
                    mean_paralogy_index=round(mean(paralogy_index),1),
                    diagnosis=ifelse(unique(orderMean)<=seg.mean$psi[2],"Single","Paralog"),
                    .groups="keep")
    
    write.table(out_meanSort_summary,paste0(outdir,"/paralogy_summary.txt"),
                quote=F,row.names=F,sep="\t")
    out_meanSort %>% group_by(qseqid) %>% arrange(paralog_percent_ignoreMissing,.by_group=TRUE) %>%
    write.table(paste0(outdir,"/paralogy_details.txt"),
                quote=F,row.names=F,sep="\t")
    
    out_meanSort_summary %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
        write.table(paste0(directory,"/paralog_list.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")
    out_meanSort_summary %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
        write.table(paste0(directory,"/singleCopy_list.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")
}

suppressMessages(suppressWarnings(require(optparse,quietly=TRUE,warn.conflicts=FALSE)))

p <- OptionParser(usage=" This script will take the CoverageStats output of VetTargets_genome.R for \n
                        many samples and detect paralogs in a targeted set of genes ")
# Add a positional argument
p <- add_option(p, c("-s","--samples"), help="<Required: list of sample names>",type="character")
p <- add_option(p, c("-d","--directory"), help="<Required: directory with output from VetTargets_genome.R>",type="character")
p <- add_option(p, c("-o","--outdir"), help="<Directory in which to write results. Defaults to current directory>",type="character",default=".")
p <- add_option(p, c("-f","--force"), help="<Force overwrite of results in outdir? Default=FALSE>",type="logical",default=FALSE)
# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(segmented,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(heatmaply,quietly=TRUE,warn.conflicts=FALSE)))

collate(samples = args$samples,
        directory = args$directory,
        outdir = args$outdir,
        force = args$force)

















