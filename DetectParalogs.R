
DetectParalogs<-function(samples,directory,outdir,force,phylogeny,ingroup){

cat("Beginning DetectParalogs with the following parameters:\n")
cat("Samples list:\t",samples,"\n")
cat("VetTargets_genome results directory:\t",directory,"\n")
cat("Output directory:\t",outdir,"\n")
cat("Force overwrite of existing results in output directory?:\t",force,"\n")
cat("input phylogeny:\t",phylogeny,"\n")
cat("ingroup samples:\t",ingroup,"\n")
cat("\n###----------------------------------------------------------------------###\n\n")


    if(dir.exists(outdir)){
        if(!force){
            stop("Specified output directory exists.\nTerminating. Use -f TRUE to overwrite.\n")
        }else{
            cat("Specified output directory exists.\n")
            cat("Will overwrite. (If it fails it's because the output file is open on your computer somewhere. Close  and rerun.)\n")
        }
    }else{
        cat("Specified output directory exists NOT. Will attempt to create it.\n")
        dir.create(outdir,recursive=TRUE)
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
    my_intercept <- min(out_meanSort$paralog_percent_ignoreMissing)
    if(my_intercept==0){
        my_formula<-formula(paralog_percent_ignoreMissing~orderMean+0)
    }else{
        my_formula<-formula(paralog_percent_ignoreMissing~orderMean)       
    }
    cat("minimum paralogy rate is",my_intercept,". Using the following formula to predict paralogs\n")
    print(my_formula)
    lm.mean<-lm(my_formula, data=out_meanSort) # Force intercept through minimum
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
        # geom_smooth(aes(x=qseqid,y=paralog_percent_ignoreMissing),method=loess) +
        scale_colour_viridis_c("Missingness (%)",option="D")+
        theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5))+
        geom_line(aes(x=orderMean,y=seg.pred))+
        geom_vline(xintercept = seg.mean$psi[2],lty=2)+
        labs(y="Paralog Rate (%)",x="Target") 
    ggsave(paste0(outdir,"/breakpoint_paralog_percent_ignoreMissing.pdf"),width=40,height=20,units = "cm")
    
    ###############
    ## HEATMAPS ###
    ###############
    
    ###############
    ### paralogy ##
    ###############
    ggplot(out_meanSort)+
        geom_raster(aes(x=qseqid,y=Sample,fill=paralog_percent_ignoreMissing))+
        scale_fill_viridis_c("Paralogy (%)",option="D")+
        labs(x="Target")+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3.5))
    ggsave(paste0(outdir,"/paralogy_heatmap.pdf"),width=28,height=20,units="cm")

    # Clustered heatmap using gplots (heatmap2()) and dendextend
    out_mat_samp_misscode<-out_meanSort %>% select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = paralog_percent_ignoreMissing, 
        values_fill = -100)
    out_mat_samp<-out_meanSort %>% select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = paralog_percent_ignoreMissing, 
        values_fill = NA)
    
    out_mat_mat<-as.data.frame(out_mat_samp_misscode)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    # order for rows
    Rowv  <- out_mat_mat %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2) %>%
        dendextend::ladderize()
    # Order for columns: We must transpose the data
    Colv  <- out_mat_mat %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>%
        set("branches_lwd", 1.2) %>%
        dendextend::ladderize()
    out_mat_mat<-as.data.frame(out_mat_samp)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    pdf(paste0(outdir,"/paralogy_heatmap2_clustered.pdf"),width=40,height=20)
    heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
              col = c(viridisLite::viridis(100)),
              Rowv = Rowv, Colv = Colv, 
              cexRow = 1 + 1/log10(nrow(out_mat_mat)),
              cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
              margins=c(10,20),
              key.title = "",
              key.xlab = "Paralogy (%)",
              key.par=list(cex.lab=2),
              adjRow = c(0,0.5),
              adjCol = c(0,0.5),
              offsetRow = 0.5,
              offsetCol = 0.5)
    dev.off()
    
    ###################################
    ##### PARALOGY with PHYLOGENY #####      
    ###################################
    if(is.null(phylogeny) | try(basename(phylogeny))=="NULL"){
        cat("no phylogeny provided. Moving on.\n")
    }else{
        cat("Phylogeny", phylogeny,"provided. Will make additional heatmap and a tanglegram.")
        ## load tree and plot
        tree<-read.tree(phylogeny)
        dend_initial<-as.dendrogram(chronos(multi2di(tree)))
        dend_initial.ord<-dend_initial %>%
            set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2)  %>%
            dendextend::ladderize(T)
        dendromatnames<-left_join(data.frame(tipnames=labels(dend_initial.ord)),
                                  data.frame(tipnames=rownames(out_mat_mat)),
                                  by="tipnames")
        # we make a matrix with rows following the order in which the samples appear in the phylogeny
        out_mat_mat_dendro_order<-out_mat_mat[match(dendromatnames$tipnames,rownames(out_mat_mat)),]
        pdf(paste0(outdir,"/paralogy_heatmap2_clustered_phylogeny.pdf"),width=40,height=20)
        heatmap.2(as.matrix(out_mat_mat_dendro_order),
                  Rowv = dend_initial.ord,
                  Colv = Colv,
                  scale = "none", col = c(viridisLite::viridis(100)),
                  trace = "none", density.info = "none",
                  cexRow = 1 + 1/log10(nrow(out_mat_mat)),
                  cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
                  margins=c(10,20),
                  key.title = "",
                  key.xlab = "Paralogy (%)",
                  key.par=list(cex.lab=2),
                  adjRow = c(0,0.5),
                  adjCol = c(0,0.5),
                  offsetRow = 0.5,
                  offsetCol = 0.5)
        dev.off()
        ##Tanglegram plotting
        obs.rf<-dist.dendlist(dendlist(Rowv,dend_initial.ord))
        pdf(paste0(outdir,"/phylogeny_vs_clustogram_paralogy_tanglegram.pdf"),width=14,height=7)
        dendextend::tanglegram(Rowv,dend_initial.ord,margin_inner = 4, lwd = 2, axes=F, 
                               main_left = "Paralogy clusters",
                               main_right= paste0("Given phylogeny\n",basename(phylogeny)),
                               highlight_distinct_edges = F,
                               common_subtrees_color_lines = F,
                               highlight_branches_lwd = F,
                               cex_main=1,
                               sub=paste0("Robinson-Foulds distance = ",obs.rf))
        dev.off()
    }

    ##################
    ### Copy Number ##
    ##################
    ggplot(out_meanSort)+
        geom_raster(aes(x=qseqid,y=Sample,fill=paralogy_index_ignoreMissing))+
        scale_fill_viridis_c("Copy number",option="D")+
        labs(x="Target")+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3.5))
    ggsave(paste0(outdir,"/Copy_number_heatmap.pdf"),width=28,height=20,units="cm")

    # Clustered heatmap using gplots (heatmap2()) and dendextend
    out_mat_samp_misscode<-out_meanSort %>% select(qseqid,Sample,paralogy_index_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = paralogy_index_ignoreMissing, 
        values_fill = -1)
    out_mat_samp<-out_meanSort %>% select(qseqid,Sample,paralogy_index_ignoreMissing) %>%
        pivot_wider(names_from = qseqid ,values_from = paralogy_index_ignoreMissing, 
        values_fill = NA)
    
    out_mat_mat<-as.data.frame(out_mat_samp_misscode)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    # order for rows
    Rowv  <- out_mat_mat %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2) %>%
        dendextend::ladderize()
    # Order for columns: We must transpose the data
    Colv  <- out_mat_mat %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>%
        set("branches_lwd", 1.2) %>%
        dendextend::ladderize()
    out_mat_mat<-as.data.frame(out_mat_samp)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    pdf(paste0(outdir,"/Copy_number_heatmap2_clustered.pdf"),width=40,height=20)
    heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
              col = c(viridisLite::viridis(100)),
              Rowv = Rowv, Colv = Colv, 
              cexRow = 1 + 1/log10(nrow(out_mat_mat)),
              cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
              margins=c(10,20),
              key.title = "",
              key.xlab = "Copy number",
              key.par=list(cex.lab=2),
              adjRow = c(0,0.5),
              adjCol = c(0,0.5),
              offsetRow = 0.5,
              offsetCol = 0.5)
    dev.off()
    
    ###############
    # missingness #
    ###############
    out_meanSort %>% mutate(qseqid=fct_reorder(qseqid,missing_percent,mean)) %>% 
        ggplot()+
        geom_raster(aes(x=qseqid,y=Sample,fill=missing_percent))+
        scale_fill_viridis_c("Missingness (%)",option="magma",direction=1)+
        labs(x="Target")+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3.5))
    ggsave(paste0(outdir,"/missingness_heatmap.pdf"),width=28,height=20,units="cm")
    # Clustered heatmap using gplots (heatmap2()) and dendextend
    out_mat_samp<-out_meanSort %>% select(qseqid,Sample,missing_percent) %>%
        pivot_wider(names_from = qseqid ,values_from = missing_percent, 
        values_fill = 100)
    
    out_mat_mat<-as.data.frame(out_mat_samp)
    rownames(out_mat_mat)<-out_mat_mat$Sample
    out_mat_mat<-out_mat_mat[,-1]
    # order for rows
    Rowv  <- out_mat_mat %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2) %>%
        dendextend::ladderize()
    # Order for columns: We must transpose the data
    Colv  <- out_mat_mat %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>%
        set("branches_lwd", 1.2) %>%
        dendextend::ladderize()

    pdf(paste0(outdir,"/missingness_heatmap2_clustered.pdf"),width=40,height=20)
    heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
              col = c(viridisLite::magma(100)),
              Rowv = Rowv, Colv = Colv, 
              cexRow = 1 + 1/log10(nrow(out_mat_mat)),
              cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
              margins=c(10,20),
              key.title = "",
              key.xlab = "Missingness (%)",
              key.par=list(cex.lab=2),
              adjRow = c(0,0.5),
              adjCol = c(0,0.5),
              offsetRow = 0.5,
              offsetCol = 0.5)
    dev.off()
    
## Per-sample plots
   #plot missingness
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/missing_percent_perSample_boxplot.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=missing_percent))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Missingness(%)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/missing_percent_perSample_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot paralogy
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_perSample_boxplot.pdf"),width=30,units = "cm")

    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Observed Paralog Rate (ignoring missing segments)",x="Sample (arranged by RSS relative to breakpoint regression)")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_perSample_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")

    # Plot residuals
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("All target loci")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_breakpoint_residuals_boxplot.pdf"),width=30,units = "cm")
    
    out_meanSort %>% mutate(Sample=fct_reorder(Sample,resid_sq,mean)) %>% 
        filter(orderMean<=seg.mean$psi[2]) %>%
        ggplot(aes(x=Sample,y=resid))+
        theme_bw()+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
        labs(y="Obs - Exp Paralog Rate for breakpoint regression with 1 inflection point")+
        ggtitle("Putative paralogs removed")
    ggsave(paste0(outdir,"/paralog_percent_ignoreMissing_breakpoint_residuals_boxplot_paralogsRemoved.pdf"),width=30,units = "cm")

    
    # Ingroup stuff #
    if(is.null(ingroup) | try(basename(ingroup))=="NULL"){
        cat("No ingroup provided. Moving on.\n")
        }else{
        cat("Using ingroup list from",ingroup,"for ingroup-specific paralog identification.\n")
        ingr<-scan(ingroup,what="character")
        out_meanSort_ingroup<-out_meanSort[out_meanSort$Sample %in% ingr,] %>%
            mutate(qseqid=fct_reorder(qseqid,paralog_percent_ignoreMissing,mean)) %>%
            arrange(paralog_percent_ignoreMissing)
        if(nrow(out_meanSort_ingroup) == 0) stop("the names in your ingroup file don't match those in the samples file")
        out_meanSort_ingroup$orderMean<-sapply(out_meanSort_ingroup$qseqid,function(x) {
            y<-which(levels(out_meanSort_ingroup$qseqid)==x)
            return(y)})
        my_intercept <- min(out_meanSort_ingroup$paralog_percent_ignoreMissing)
        if(my_intercept==0){
            my_formula<-formula(paralog_percent_ignoreMissing~orderMean+0)
        }else{
            my_formula<-formula(paralog_percent_ignoreMissing~orderMean)       
        }
        cat("INGROUP minimum paralogy rate is",my_intercept,". Using the following formula to predict paralogs\n")
        print(my_formula)

        lm.mean.ingroup<-lm(my_formula, data=out_meanSort_ingroup) # Force intercept through minimum
        seg.mean.ingroup<-segmented(lm.mean.ingroup,npsi=1)
        out_meanSort_ingroup$seg.pred<-predict(seg.mean.ingroup,newdata=data.frame(orderMean=out_meanSort_ingroup$orderMean))
        out_meanSort_ingroup$resid<-out_meanSort_ingroup$paralog_percent_ignoreMissing - out_meanSort_ingroup$seg.pred
        out_meanSort_ingroup$resid_sq<-out_meanSort_ingroup$resid^2      
    
        ##plot breakpoint analysis
        out_meanSort_ingroup %>%
            ggplot(aes(x=qseqid,y=paralog_percent_ignoreMissing))+
            theme_bw()+
            geom_point(aes(colour=missing_percent),alpha=0.5)+
            # geom_smooth(aes(x=qseqid,y=paralog_percent_ignoreMissing),method=loess) +
            scale_colour_viridis_c("Missingness (%)",option="D")+
            theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5))+
            geom_line(aes(x=orderMean,y=seg.pred))+
            geom_vline(xintercept = seg.mean.ingroup$psi[2],lty=2)+
            labs(y="Paralog Rate (%)",x="Target") 
        ggsave(paste0(outdir,"/breakpoint_paralog_percent_ignoreMissing_ingroup.pdf"),width=40,height=20,units = "cm")

        out_meanSort_summary_ingroup<-out_meanSort_ingroup %>% group_by(qseqid) %>%
        summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                    mean_missing_percent=round(mean(missing_percent),1),
                    mean_paralogy_index=round(mean(paralogy_index),2),
                    diagnosis=ifelse(unique(orderMean)<=seg.mean.ingroup$psi[2],"Single","Paralog"),
                    .groups="keep")
        write.table(out_meanSort_summary_ingroup,paste0(outdir,"/DetectParalogs_results_summarised_ingroup.txt"),
                    quote=F,row.names=F,sep="\t")
        cat("Using ingroup samples -- Identified",sum(out_meanSort_summary_ingroup$diagnosis=="Paralog"),"paralogs!\n")

        out_meanSort_summary_ingroup %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
        write.table(paste0(outdir,"/paralog_list_ingroup.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")
        out_meanSort_summary_ingroup %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
        write.table(paste0(outdir,"/singleCopy_list_ingroup.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")
    }
    # ----
    
    #####################
    ### DONE PLOTTING ###
    #####################
    
    out_meanSort %>% group_by(qseqid) %>% arrange(paralog_percent_ignoreMissing,.by_group=TRUE) %>% as.data.frame() %>%
        write.table(paste0(outdir,"/DetectParalogs_results_full.txt"),quote=F,row.names=F,sep="\t")

    out_meanSort_summary<-out_meanSort %>% group_by(qseqid) %>%
        summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                    mean_missing_percent=round(mean(missing_percent),1),
                    mean_paralogy_index=round(mean(paralogy_index),2),
                    mean_paralogy_index_ignoreMissing=round(mean(paralogy_index_ignoreMissing),2),
                    diagnosis=ifelse(unique(orderMean)<=seg.mean$psi[2],"Single","Paralog"),
                    .groups="keep")
    write.table(out_meanSort_summary,paste0(outdir,"/DetectParalogs_results_summarised.txt"),
                quote=F,row.names=F,sep="\t")
    cat("Using all samples -- Identified",sum(out_meanSort_summary$diagnosis=="Paralog"),"paralogs!\n")
    
    out_meanSort_summary %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
        write.table(paste0(outdir,"/paralog_list.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")
    out_meanSort_summary %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
        write.table(paste0(outdir,"/singleCopy_list.txt"),
                    quote=F,row.names=F,col.names=F,sep="\t")

}



suppressMessages(suppressWarnings(require(optparse,quietly=TRUE,warn.conflicts=FALSE)))

p <- OptionParser(usage=" This script will take the CoverageStats output of VetTargets_genome.R for \n
                        many samples and detect paralogs in a targeted set of genes ")
# Add a positional argument
p <- add_option(p, c("-s","--samples"), help="<Required: list of sample names>",type="character")
p <- add_option(p, c("-d","--directory"), help="<Required: directory with output from VetTargets_genome.R>",type="character",default=paste0(getwd(),"/DetectParalogs_results"))
p <- add_option(p, c("-o","--outdir"), help="<Required: directory in which to write results>",type="character")
p <- add_option(p, c("-f","--force"), help="<Force overwrite of results in outdir? Default=FALSE>",type="logical",default=FALSE)
p <- add_option(p, c("-p","--phylogeny"), help="<Rooted tree in Newick format. If provided, will make an additional paralogy heatmap with this tree instead of the cluster dendrogram. All tip labels need to match those in the samples file.>",type="character",default=NULL)
p <- add_option(p, c("-i","--ingroup"), help="<File listing 'ingroup' samples. This is useful if you have several outgroup taxa, which often have different paralogy patterns (especially if they were used to design the target set). \nA separate paralog detection analysis will be conducted using only the ingroup samples.>",type="character",default=NULL)
# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(segmented,quietly=TRUE,warn.conflicts=FALSE)))
# suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dplyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(forcats,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(gplots,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dendextend,quietly=TRUE,warn.conflicts=FALSE)))

if(!is.null(args$phylogeny) & args$phylogeny != "NULL") suppressMessages(suppressWarnings(require(ape,quietly=TRUE,warn.conflicts=FALSE)))


try(DetectParalogs(samples = args$samples,
        directory = args$directory,
        outdir = args$outdir,
        force = args$force,
        phylogeny = args$phylogeny,
        ingroup = args$ingroup))

















