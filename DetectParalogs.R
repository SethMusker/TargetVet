
DetectParalogs<-function(samples,directory,outdir,force,phylogeny=NULL,ingroup=NULL,residual_cutoff){
  
  cat("Beginning DetectParalogs with the following parameters:\n")
  cat("Samples list:\t",samples,"\n")
  cat("VetTargets_genome results directory:\t",directory,"\n")
  cat("Output directory:\t",outdir,"\n")
  cat("Force overwrite of existing results in output directory?:\t",force,"\n")
  cat("input phylogeny:\t",phylogeny,"\n")
  cat("ingroup samples:\t",ingroup,"\n")
  cat("\n###-------------------------------------------------###\n\n")
  
  
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
  
  samples<-suppressMessages(scan(samples,what="character"))
  
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
    mutate(qseqid=fct_reorder(qseqid,paralog_percent_ignoreMissing,mean)) 
  # %>% arrange(paralog_percent_ignoreMissing)
  
  #get the order and do segmented
  out_meanSort$orderMean<-sapply(out_meanSort$qseqid,function(x) {
    y<-which(levels(out_meanSort$qseqid)==x)
    return(y)})
  
  out_meanSort<-out_meanSort %>% 
    arrange(orderMean)
  
  ## do piecewise constant aka step-function model (package 'cumSeg')
  cat("Fitting piecewise constant a.k.a step-function model with # breakpoints = 1.\n")
  step_mod<-jumpoints(out_meanSort$paralog_percent_ignoreMissing,
                      out_meanSort$orderMean,k=1,output = "2")
  ## TODO: fix this! Predicted values don't go in right.
  out_meanSort$step.k1.pred<-step_mod$fitted.values
  step.psi<-step_mod$psi
  cat("Breakpoint estimate =",step.psi,"meaning the step-function estimated number of paralogs is\n\t",
      max(out_meanSort$orderMean),"-",step.psi,"=",max(out_meanSort$orderMean)-step.psi,".\n\n")
  
  
  ## do n-parameter logistic regression with automated selection of best n  
  cat("Fitting n-parameter logistic regression (package 'nplr') with automated selection of best n.\n")    
  nplr_mod<-nplr(x=out_meanSort$orderMean,
                 y=out_meanSort$paralog_percent_ignoreMissing/100,
                 npars = "all",
                 useLog = F)
  cat("The best-fitting model had",getPar(nplr_mod)$npar,"parameters:\n")
  print(getPar(nplr_mod)$params)
  cat("\n\n")
  # nplr_mod.psi<-getInflexion(nplr_mod)[1][,1]
  # cat("Inflexion point estimate =",nplr_mod.psi,"meaning the n-parameter logistic regression estimated number of paralogs is\n\t",
  #     max(out_meanSort$orderMean),"-",nplr_mod.psi,"=",max(out_meanSort$orderMean)-nplr_mod.psi,".\n")
  
  ## Rather define psi as NP50 (x at which line crosses y=50%)
  nplr_mod.psi<-getEstimates(nplr_mod,0.5)$x
  cat("NP50 estimate =",round(nplr_mod.psi,0),"meaning the n-parameter logistic regression estimated number of paralogs is\n\t",
      max(out_meanSort$orderMean),"-",round(nplr_mod.psi,0),"=",max(out_meanSort$orderMean)-round(nplr_mod.psi,0),".\n")
  
  ## find anomalous samples based on how well they conform to the overall paralogy patterns
  ## get the predicted values of the overall model
  predict_nplr<-function(model,X){
    pars<-getPar(model)$params
    bottom <- pars[1][1,]
    top <- pars[2][1,]
    xmid <- pars[3][1,]
    scal <- pars[4][1,]
    s <- pars[5][1,]
    pred<-bottom + (top - bottom)/(1 + 10^((xmid - X) * scal))^s
    return(pred)
  }
  out_meanSort$nplr_mod.pred<-predict_nplr(nplr_mod,out_meanSort$orderMean)*100
  out_meanSort$nplr_mod.pred_resid<-out_meanSort$paralog_percent_ignoreMissing - out_meanSort$nplr_mod.pred 
  ## We expect the residuals to equal zero, so we force the intercept to zero
  resid_model<-lm(nplr_mod.pred_resid~Sample+0,data=out_meanSort)
  resid_model_coefs<-as.data.frame(summary(resid_model)$coefficients)
  names(resid_model_coefs)<-c("Mean_residual_paralogy","SE","t_value","p_value")
  resid_model_coefs$Sample<-factor(gsub("Sample","",row.names(resid_model_coefs)))
  residual_deviants<-resid_model_coefs[abs(resid_model_coefs$Mean_residual_paralogy)>residual_cutoff,"Sample"]
  
  
  ## Make some objects for plotting
  ## Get predicted values and 95% CIs
  newx <- getXcurve(nplr_mod)
  newy <- getYcurve(nplr_mod)
  bounds <- nplr:::.confInt(getStdErr(nplr_mod)["stdErr"], getY(nplr_mod), newy)
  nplr_mod_estimates<-data.frame(x=newx,y=newy*100,y.975=bounds$hi*100,y.025=bounds$lo*100)
  nplr_mod_estimates$y.975[nplr_mod_estimates$y.975>100]<-NA
  nplr_mod_estimates$y.025[nplr_mod_estimates$y.025<0]<-NA
  
  # inflexion_estimates<-getEstimates(nplr_mod,getInflexion(nplr_mod)[2][,1])
  inflexion_estimates<-getEstimates(nplr_mod,0.5)
  
  ## fit separate nplr models to each sample, check for concordance of inflexion estimates
  get_np_fit<-function(x,y){
    mod_temp<-suppressMessages(nplr(x=x,
                                    y=y/100,
                                    npars = "all",
                                    useLog = F))
    newx <- getXcurve(mod_temp)
    newy <- getYcurve(mod_temp)
    para_quarts<-getEstimates(mod_temp,c(0.25,0.5,0.75))[,c(1,3)]
    if(any(para_quarts$y != c(0.25,0.5,0.75))){
      cat("Curve does not pass through some quartiles. Setting estimates to NA.\n")
      para_quarts$y<-c(0.25,0.5,0.75)
      para_quarts$x<-NA
    }
    
    return(list(curves=data.frame(x=newx,y=newy*100),
                model=mod_temp,
                inflexions=getInflexion(mod_temp),
                paralogy_quartiles=para_quarts))
  }
  
  out_meanSort_sampleSplit<-split(out_meanSort,~Sample)
  for(i in names(out_meanSort_sampleSplit)){
    out_meanSort_sampleSplit[[i]] <- out_meanSort_sampleSplit[[i]] %>% 
      mutate(qseqid=fct_reorder(qseqid,paralog_percent_ignoreMissing,mean)) %>%
      arrange(paralog_percent_ignoreMissing)
  }
  cat("Now fitting nplr models separately to each sample.\n")
  myfits<-lapply(out_meanSort_sampleSplit,function(S) get_np_fit(x=S$orderMean,y=S$paralog_percent_ignoreMissing))
  mycurves<-lapply(myfits,function(x) return(x$curves))
  myinflexions<-lapply(myfits,function(x) return(x$inflexions[1][,1]))
  myinflexions_y<-lapply(myfits,function(x) return(x$inflexions[2][,1]*100))
  myNP50<-lapply(myfits,function(x) return(x$paralogy_quartiles[2,"x"]))
  
  # TODO: use quartiles, esp. median i.e. where line crosses 50% paralogy, instead of inflexion points
  # how to deal with flat curves? --  for now, set estimates to NA. We anyway do the residual analysis to find outliers.
  
  for(i in names(mycurves)){
    mycurves[[i]]$Sample<-i
  }
  mycurves_combined<-do.call(rbind,mycurves)
  row.names(mycurves_combined)<-1:nrow(mycurves_combined)
  mycurves_combined$y[mycurves_combined$y<0]<-0
  mycurves_combined$y[mycurves_combined$y>100]<-100
  
  inflexions_df<-data.frame(Sample=names(myinflexions),
                            inflexion=do.call(c,myinflexions),
                            inflexion_y=do.call(c,myinflexions_y),
                            NP50=do.call(c,myNP50))
  
  ## Get deviation from 'normality' as per default of qqplot (only if there are >= 5 samples)
  # cat("Now looking for anomalous (i.e. outlier) samples based on their estimated number of single-copy targets.\n")
  # if(nrow(inflexions_df)>=5){
  #   cat("Number of samples is >= 5. Will look for anomalous samples based on 'robust normality' test.\n")
  #   myqq<-qqnorm(inflexions_df$inflexion,plot.it = F)
  #   ## the line drawn by qqline goes through the 1st and 3rd quartiles (more robust than sigma*x + mean)
  #   ydiff<-dist(quantile(inflexions_df$inflexion,c(1/4,3/4)))
  #   xdiff<-dist(quantile(myqq$x,c(1/4,3/4)))
  #   qq_slope<-ydiff/xdiff
  #   qq_mid<-mean(quantile(inflexions_df$inflexion,c(1/4,3/4)))
  #   ## y = mx + c
  #   ## inflexion(exp) = qq_slope*myqq$x + qq_mid
  #   # TODO: output qqplot
  #   # qqnorm(inflexions_df$inflexion);qqline(inflexions_df$inflexion, col = 2)
  #   # lines(myqq$x,qq_slope*myqq$x + qq_mid)
  #   expected<-qq_slope*myqq$x + qq_mid
  #   observed<-myqq$y
  #   inflexions_df$qq_deviation<-observed-expected
  #   # Flag as an outlier if it (absolutely) deviates from the expectation by more than 10% of the total number of genes
  #   inflexions_df$big_deviant<-abs(inflexions_df$qq_deviation)/(max(out_meanSort$orderMean)/10) > 1
  #   if(any(inflexions_df$big_deviant)){
  #     cat("WARNING: some samples seem to have anomalous paralogy patterns!\n")
  #     cat("Specifically, the 'robust' absolute deviation from normality of their estimated number of single-copy targets exceeds 10% of the total number of genes in the data set.\n")
  #     write.table(data.frame(x=inflexions_df$sample[which(inflexions_df$big_deviant)]),
  #                 file="Anomalous_samples.csv",quote=F,row.names=F,col.names=F)
  #     cat("The names of these samples are written to the file 'Anomalous_samples.csv'.\n")
  #   }
  #   
  #   
  # } else{
  #   # if too few samples, simply identify big deviants by whether their inflexions fall outside possible bounds
  #   inflexions_df$big_deviant<-inflexions_df$inflexion > max(out_meanSort$orderMean) | inflexions_df$inflexion < 1
  #   if(any(inflexions_df$big_deviant)){
  #     cat("WARNING: some samples seem to have anomalous paralogy patterns!\n")
  #     cat("Specifically, their estimated number of single-copy targets either exceeds the total number of genes in the data set, or is effectively zero.\n")
  #     write.table(data.frame(x=inflexions_df$sample[which(inflexions_df$big_deviant)]),
  #                 file="Anomalous_samples.csv",quote=F,row.names=F,col.names=F)
  #     cat("The names of these samples are written to the file 'Anomalous_samples.csv'.\n")
  #   }
  # }
  # much simpler than above: define deviants based on their mean residuals
  inflexions_df$big_deviant<-inflexions_df$Sample %in% residual_deviants
  if(length(residual_deviants)>0){
    cat("WARNING:",length(residual_deviants),"samples seem to have anomalous paralogy patterns!\n")
    cat("Specifically, their absolute mean observed minus expected paralogy rate exceeds",residual_cutoff,"%.\n")
    write.table(data.frame(x=inflexions_df$Sample[which(inflexions_df$big_deviant)]),
                file=paste0(outdir,"/Paralogy_anomalous_samples.txt"),quote=F,row.names=F,col.names=F)
    cat("The names of these samples are written to the file 'Paralogy_anomalous_samples.txt'.\n")
  }
    # make data for samplewise plot
  inflexions_df$plot_point_x<-ifelse(inflexions_df$NP50 > max(out_meanSort$orderMean),
                                     yes=max(out_meanSort$orderMean),
                                     no=ifelse(inflexions_df$NP50 < 1,1,inflexions_df$NP50)) # plot points at their NP50
  inflexions_df$plot_point_x[is.na(inflexions_df$plot_point_x)]<-max(out_meanSort$orderMean)/2
  inflexions_df$plot_point_y<-NA
  for(i in inflexions_df$Sample){
    inflexions_df[inflexions_df$Sample==i,]$plot_point_y<-ifelse(inflexions_df[inflexions_df$Sample==i,]$inflexion_y <= 100 & inflexions_df[inflexions_df$Sample==i,]$inflexion_y > 0,
                                                                 yes=inflexions_df[inflexions_df$Sample==i,]$inflexion_y,
                                                                 no=mean(out_meanSort[out_meanSort$Sample==i,]$paralog_percent_ignoreMissing))  
  }                                
  
  inflexions_df$label<-ifelse(inflexions_df$big_deviant,inflexions_df$Sample,NA)
  inflexions_df$label_hjust<-ifelse(inflexions_df$inflexion > max(out_meanSort$orderMean),
                                    yes=1,
                                    no=ifelse(inflexions_df$inflexion < 1,0,
                                              ifelse(inflexions_df$inflexion < nplr_mod.psi,1,0)))
  
  # Write out results
  inflexions_df_out<-inflexions_df
  inflexions_df_out$inflexion<-round(inflexions_df_out$inflexion,1)
  inflexions_df_out<-left_join(inflexions_df_out,resid_model_coefs,"Sample")
  write.table(inflexions_df_out,
              paste0(outdir,"/Paralogy_estimates_samplewise.csv"),
              quote=F,row.names=F,sep=",")
  cat("The results of the sample-wise paralogy and anomalous sample detection analysis are written to 'Paralogy_estimates_samplewise.csv'.\n")
  

  mycurves_combined<-left_join(mycurves_combined,
                               inflexions_df[,c("Sample","big_deviant")],
                               "Sample")
  
  
  #####################
  ### BEGIN PLOTTING ##
  #####################
  # colours from https://davidmathlogic.com/colorblind/#%231A85FF-%23D41159
  step_nplr_plot<-ggplot()+
    theme_bw()+
    labs(y="Paralog Rate (%)",x="Target") +
    theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5,hjust=1),
          panel.grid =element_blank())+
    ggtitle("n-parameter logistic regression (magenta) and step-function (blue) model predictions.")+
    geom_bin2d(aes(x=qseqid,
                   y=paralog_percent_ignoreMissing),
               show.legend = F,data=out_meanSort)+
    scale_fill_gradient(low="gray90",high="black")+
    ylim(-5,105)+
    geom_line(aes(x=x,y=y),data=nplr_mod_estimates,size=1.5,col="#D41159")+ # nplr fit
    geom_line(aes(x=x,y=y.025),data=nplr_mod_estimates,colour="#D41159",lty=3,size=1)+ # lower CI
    geom_line(aes(x=x,y=y.975),data=nplr_mod_estimates,colour="#D41159",lty=3,size=1)+ # upper CI
    # geom_vline(xintercept = inflexion_estimates$x,lty=2,size=1,col="#D41159")+ 
    # geom_vline(xintercept = inflexion_estimates$x.025,lty=2,colour="#D41159",size=1)+
    # geom_vline(xintercept = inflexion_estimates$x.975,lty=2,colour="#D41159",size=1)+
    geom_segment(aes(x = x,xend=x,y=0,yend=100),
                 data=inflexion_estimates,
                 lty=2,
                 size=1,
                 col="#D41159")+ 
    geom_segment(aes(x = x.025,xend=x.025,y=0,yend=100),
                 data=inflexion_estimates,
                 lty=2,
                 size=0.75,
                 col="#D41159")+ 
    geom_segment(aes(x = x.975,xend=x.975,y=0,yend=100),
                 data=inflexion_estimates,
                 lty=2,
                 size=0.75,
                 col="#D41159")+ 
    geom_line(aes(x=orderMean,y=step.k1.pred),data=out_meanSort,colour="#1A85FF",size=1.5)+ # step-model fit
    geom_segment(aes(x = x,xend=x,y=0,yend=100),
                 data=data.frame(x=step.psi),
                 lty=2,
                 size=1,
                 col="#1A85FF") +
    geom_hline(yintercept=50,lty=2,size=0.5)   
  # geom_vline(xintercept = step.psi,lty=2,colour="#1A85FF",size=1)      
  ggsave(paste0(outdir,"/Paralog_detection_plot.pdf"),
         plot=step_nplr_plot,
         width=40,height=20,units = "cm")
  
  samplewise_nplr_plot<-ggplot()+
    theme_bw()+
    labs(y="Paralog Rate (%)",x="Target") +
    theme(axis.text.x = element_text(angle = 90,size=3.5,vjust=0.5,hjust=1),
          panel.grid =element_blank())+
    ggtitle("Sample-wise n-parameter logistic regression model predictions.")+
    geom_bin2d(aes(x=qseqid,
                   y=paralog_percent_ignoreMissing),
               show.legend = F,data=out_meanSort)+
    scale_fill_gradient(low="gray90",high="black")+
    ylim(-5,105)+
    geom_line(aes(x=x,y=y),data=nplr_mod_estimates,size=1,col="#D41159")+
    # add the full model estimates (inflexion and curve)
    geom_segment(aes(x = x,xend=x,y=0,yend=100),
                 data=inflexion_estimates,
                 lty=2,
                 size=1,
                 col="#D41159")+ 
    geom_line(aes(x=x,y=y,group=Sample,colour=big_deviant,alpha=big_deviant),
              data=mycurves_combined,
              show.legend = F)+
    # add samplewise curves and inflexion points
    # geom_segment(aes(x=plot_point_x,xend=plot_point_x,
    #                  y=0,yend=inflexion_y,colour=big_deviant,
    #                  alpha=big_deviant),
    #              data=inflexions_df,
    #              lty=2,
    #              show.legend = F)+
    scale_alpha_manual(values=c(0.5,1))+
    # label outliers
    geom_text(aes(x=plot_point_x,y=plot_point_y,label=label,hjust=label_hjust),
              data=inflexions_df,na.rm=T)+
    # geom_point(aes(x=plot_point_x,y=inflexion_y,colour=big_deviant),
    #            data=inflexions_df,show.legend = F)+
    scale_colour_manual(values=c("blue","#FFC107"))+
    geom_hline(yintercept=50,lty=2,size=0.5) 
  ggsave(paste0(outdir,"/Paralog_detection_plot_samplewise.pdf"),
         plot=samplewise_nplr_plot,
         width=40,height=20,units = "cm")
  
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
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5))
  ggsave(paste0(outdir,"/Paralogy_heatmap.pdf"),width=28,height=20,units="cm")
  
  # Clustered heatmap using gplots (heatmap2()) and dendextend
  out_mat_samp_misscode<-out_meanSort %>% 
    select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
    pivot_wider(names_from = qseqid ,values_from = paralog_percent_ignoreMissing, 
                values_fill = -100)
  out_mat_samp<-out_meanSort %>% 
    select(qseqid,Sample,paralog_percent_ignoreMissing) %>%
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
  pdf(paste0(outdir,"/Paralogy_heatmap_clustered.pdf"),width=40,height=20)
  heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
            col = c(viridisLite::viridis(100)),
            Rowv = Rowv, Colv = Colv, 
            cexRow = 1 + 1/log10(nrow(out_mat_mat)),
            cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
            keysize = 0.75,
            margins=c(10,30),
            key.title = "",
            key.xlab = "Paralogy (%)",
            key.par=list(cex.lab=3),
            adjRow = c(0,0.5),
            adjCol = c(1,0.5), ## TODO: Need to fix these!
            offsetRow = 0,
            offsetCol = 0)
  dev.off()
  
  ###################################
  ##### PARALOGY with PHYLOGENY #####      
  ###################################
  # if(is.null(phylogeny) | try(basename(phylogeny))=="NULL"){
  if(is.null(phylogeny)){
    cat("\n")
    # cat("No phylogeny provided. Moving on.\n")
  }else{
    cat("\nPhylogeny", phylogeny,"provided. Will make additional heatmap and a tanglegram.\n")
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
    pdf(paste0(outdir,"/Paralogy_heatmap_clustered_phylogeny.pdf"),width=40,height=20)
    heatmap.2(as.matrix(out_mat_mat_dendro_order),
              Rowv = dend_initial.ord,
              Colv = Colv,
              scale = "none", col = c(viridisLite::viridis(100)),
              trace = "none", density.info = "none",
              cexRow = 1 + 1/log10(nrow(out_mat_mat)),
              cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
              keysize = 0.75,
              margins=c(10,30),
              key.title = "",
              key.xlab = "Paralogy (%)",
              key.par=list(cex.lab=2),
              adjRow = c(0,0.5),
              adjCol = c(0,0.5),
              offsetRow = 0,
              offsetCol = 0)
    dev.off()
    ##Tanglegram plotting
    obs.rf<-dist.dendlist(dendlist(Rowv,dend_initial.ord))
    pdf(paste0(outdir,"/Phylogeny_vs_clustogram_paralogy_tanglegram.pdf"),width=14,height=7)
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
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5))
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
  # Rowv  <- out_mat_mat %>% dist %>% hclust %>% as.dendrogram %>%
  #   set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2) %>%
  #   dendextend::ladderize()
  # # Order for columns: We must transpose the data
  # Colv  <- out_mat_mat %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  #   set("branches_k_color", k = 4) %>%
  #   set("branches_lwd", 1.2) %>%
  #   dendextend::ladderize()
  out_mat_mat<-as.data.frame(out_mat_samp)
  rownames(out_mat_mat)<-out_mat_mat$Sample
  out_mat_mat<-out_mat_mat[,-1]
  pdf(paste0(outdir,"/Copy_number_heatmap_clustered.pdf"),width=40,height=20)
  heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
            col = c(viridisLite::viridis(100)),
            Rowv = Rowv, Colv = Colv, 
            cexRow = 1 + 1/log10(nrow(out_mat_mat)),
            cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
            keysize = 0.75,
            margins=c(10,30),
            key.title = "",
            key.xlab = "Copy number",
            key.par=list(cex.lab=2),
            adjRow = c(0,0.5),
            adjCol = c(1,0.5),
            offsetRow = 0,
            offsetCol = 0,
            main="Note: dendrograms are calculated from the paralogy matrix.")
  dev.off()
  
  ###############
  # missingness #
  ###############
  out_meanSort %>% mutate(qseqid=fct_reorder(qseqid,missing_percent,mean)) %>% 
    ggplot()+
    geom_raster(aes(x=qseqid,y=Sample,fill=missing_percent))+
    scale_fill_viridis_c("Missingness (%)",option="magma",direction=1)+
    labs(x="Target")+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5))
  ggsave(paste0(outdir,"/Missingness_heatmap.pdf"),width=28,height=20,units="cm")
  # Clustered heatmap using gplots (heatmap2()) and dendextend
  out_mat_samp<-out_meanSort %>% select(qseqid,Sample,missing_percent) %>%
    pivot_wider(names_from = qseqid ,values_from = missing_percent, 
                values_fill = 100)
  
  out_mat_mat<-as.data.frame(out_mat_samp)
  rownames(out_mat_mat)<-out_mat_mat$Sample
  out_mat_mat<-out_mat_mat[,-1]
  # order for rows
  # Rowv  <- out_mat_mat %>% dist %>% hclust %>% as.dendrogram %>%
  #   set("branches_k_color", k = 4) %>% set("branches_lwd", 1.2) %>%
  #   dendextend::ladderize()
  # # Order for columns: We must transpose the data
  # Colv  <- out_mat_mat %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  #   set("branches_k_color", k = 4) %>%
  #   set("branches_lwd", 1.2) %>%
  #   dendextend::ladderize()
  
  pdf(paste0(outdir,"/Missingness_heatmap_clustered.pdf"),width=40,height=20)
  heatmap.2(as.matrix(out_mat_mat), scale = "none", trace = "none", density.info = "none",
            col = c(viridisLite::magma(100)),
            Rowv = Rowv, Colv = Colv, 
            cexRow = 1 + 1/log10(nrow(out_mat_mat)),
            cexCol = 0.1+1/sqrt(nrow(out_mat_mat)),
            keysize = 0.75,
            margins=c(10,30),
            key.title = "",
            key.xlab = "Missingness (%)",
            key.par=list(cex.lab=2),
            adjRow = c(0,0.5),
            adjCol = c(1,0.5),
            offsetRow = 0,
            offsetCol = 0,
            main="Note: dendrograms are calculated from the paralogy matrix.")
  dev.off()
  
  ## Per-sample plots
  #plot missingness
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    ggplot(aes(x=Sample,y=missing_percent))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
    labs(y="Missingness(%)",x="Sample")+
    ggtitle("All targets.")
  ggsave(paste0(outdir,"/Missingness_boxplot_samplewise.pdf"),width=30,units = "cm")
  
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    filter(orderMean<=nplr_mod.psi) %>%
    ggplot(aes(x=Sample,y=missing_percent))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
    labs(y="Missingness(%)",x="Sample")+
    ggtitle("Putative paralogs removed.")
  ggsave(paste0(outdir,"/Missingness_boxplot_samplewise_paralogsRemoved.pdf"),width=30,units = "cm")
  
  # Plot paralogy
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
    labs(y="Paralog Rate (%)",x="Sample")+
    ggtitle("All targets.")
  ggsave(paste0(outdir,"/Paralogy_boxplot_samplewise.pdf"),width=30,units = "cm")
  
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    filter(orderMean<=nplr_mod.psi) %>%
    ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1))+
    labs(y="Paralog Rate (%)",x="Sample")+
    ggtitle("Putative paralogs removed.")
  ggsave(paste0(outdir,"/Paralogy_boxplot_samplewise_paralogsRemoved.pdf"),width=30,units = "cm")
  
  # # Plot residuals
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>%
    ggplot(aes(x=Sample,y=nplr_mod.pred_resid))+
    # geom_violin(scale="width")+  
    geom_boxplot(outlier.size = 0.5,outlier.shape = 1)+  
    geom_hline(yintercept=0,lty=2,size=1)+
    geom_hline(yintercept=residual_cutoff,lty=2,size=0.5)+
    geom_hline(yintercept=-residual_cutoff,lty=2,size=0.5)+
    geom_point(aes(x=Sample,y=Mean_residual_paralogy),data=resid_model_coefs,colour="red")+
    theme_bw()+
    labs(y="Obs - Exp paralogy (%)",x="Sample") +
    theme(axis.text.x = element_text(angle = 90,size=5.5,vjust=0.5,hjust=1),
          panel.grid =element_blank())
  ggsave(paste0(outdir,"/Paralogy_residuals_boxplot_samplewise.pdf"),width=30,units = "cm")
  
  
  # Ingroup stuff #
  if(is.null(ingroup)){
    cat("No ingroup provided. Moving on.\n")
  }else{
    cat("Using ingroup list from",ingroup,"for ingroup-specific paralog identification.\n")
    ingr<-suppressMessages(scan(ingroup,what="character"))
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
    write.table(out_meanSort_summary_ingroup,paste0(outdir,"/DetectParalogs_results_summarised_ingroup.csv"),
                quote=F,row.names=F,sep=",")
    cat("Using ingroup samples -- Identified",sum(out_meanSort_summary_ingroup$diagnosis=="Paralog"),"paralogs!\n")
    
    out_meanSort_summary_ingroup %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
      write.table(paste0(outdir,"/Paralog_list_ingroup.csv"),
                  quote=F,row.names=F,col.names=F,sep=",")
    out_meanSort_summary_ingroup %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
      write.table(paste0(outdir,"/SingleCopy_list_ingroup.csv"),
                  quote=F,row.names=F,col.names=F,sep=",")
  }
  # ----
  
  #####################
  ### DONE PLOTTING ###
  #####################
  
  out_meanSort %>% group_by(qseqid) %>% 
    arrange(paralog_percent_ignoreMissing,.by_group=TRUE) %>% 
    as.data.frame() %>%
    write.table(paste0(outdir,"/DetectParalogs_results_full.csv"),
                quote=F,row.names=F,sep=",")
  
  out_meanSort_summary<-out_meanSort %>% group_by(qseqid) %>%
    summarise(  mean_paralog_percent_ignoreMissing=round(mean(paralog_percent_ignoreMissing),1),
                mean_missing_percent=round(mean(missing_percent),1),
                mean_paralogy_index=round(mean(paralogy_index),2),
                mean_paralogy_index_ignoreMissing=round(mean(paralogy_index_ignoreMissing),2),
                diagnosis=ifelse(unique(orderMean)<=nplr_mod.psi,"Single","Paralog"),
                .groups="keep")
  write.table(out_meanSort_summary,paste0(outdir,"/DetectParalogs_results_summarised.csv"),
              quote=F,row.names=F,sep=",")
  cat("Using all samples -- Identified",sum(out_meanSort_summary$diagnosis=="Paralog"),"paralogs!\n")
  
  out_meanSort_summary %>% filter(diagnosis=="Paralog") %>% select(qseqid) %>%
    write.table(paste0(outdir,"/Paralog_list.csv"),
                quote=F,row.names=F,col.names=F,sep=",")
  out_meanSort_summary %>% filter(diagnosis=="Single") %>% select(qseqid) %>%
    write.table(paste0(outdir,"/SingleCopy_list.csv"),
                quote=F,row.names=F,col.names=F,sep=",")
  
  cat("\nDetectParalogs.R finished!\n")
  
}



suppressMessages(suppressWarnings(require(optparse,quietly=TRUE,warn.conflicts=FALSE)))

p <- OptionParser(usage=" This script will take the CoverageStats output of VetTargets_genome.R for \n
                        many samples and detect paralogs in a targeted set of genes. ")
# Add a positional arguments
p <- add_option(p, c("-s","--samples"), help="<Required: file listing sample names, one per line>",type="character")
p <- add_option(p, c("-d","--directory"), help="<Required: directory with output from VetTargets_genome.R>",type="character")
p <- add_option(p, c("-o","--outdir"), help="<Directory in which to write results. Defaults to ./DetectParalogs_results>",type="character",default=paste0(getwd(),"/DetectParalogs_results"))
p <- add_option(p, c("-f","--force"), help="<Force overwrite of results in outdir? Default=FALSE>",type="logical",default=FALSE)
p <- add_option(p, c("-p","--phylogeny"), help="<Rooted tree in Newick format. If provided, will make an additional paralogy heatmap with this tree instead of the cluster dendrogram. All tip labels need to match those in the samples file.>",type="character",default=NULL)
p <- add_option(p, c("-i","--ingroup"), help="<File listing 'ingroup' samples. This is useful if you have several outgroup taxa, which often have different paralogy patterns (especially if they were used to design the target set). \nA separate paralog detection analysis will be conducted using only the ingroup samples.>",type="character",default=NULL)
p <- add_option(p, c("-r","--residual_cutoff"), help="<How to identify anomalous samples? If their absolute mean observed minus expected paralogy exceeds this value. See documentation for details.>",type="numeric",default=10)
# parse
args<-parse_args(p)

if(any(c(is.null(args$samples),is.null(args$directory)))) {
  print_help(p)
  stop("Some required arguments are missing.")
}

# suppressMessages(suppressWarnings(require(segmented,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(cumSeg,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(nplr,quietly=TRUE,warn.conflicts=FALSE)))
# suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(tidyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dplyr,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(forcats,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(gplots,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(dendextend,quietly=TRUE,warn.conflicts=FALSE)))

# if(!is.null(args$phylogeny) & args$phylogeny != "NULL") suppressMessages(suppressWarnings(require(ape,quietly=TRUE,warn.conflicts=FALSE)))
if(!is.null(args$phylogeny)) suppressMessages(suppressWarnings(require(ape,quietly=TRUE,warn.conflicts=FALSE)))


try(DetectParalogs(samples = args$samples,
                   directory = args$directory,
                   outdir = args$outdir,
                   force = args$force,
                   phylogeny = args$phylogeny,
                   ingroup = args$ingroup,
                   residual_cutoff=args$residual_cutoff))

sink(paste0(args$outdir,"/Warnings.txt"))
warnings()
sink()
cat("Warnings are written to warnings.txt. These are usually harmless.\n")
















