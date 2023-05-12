plot_missingness_samplewise<-function(){
out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    ggplot(aes(x=Sample,y=missing_percent))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
    labs(y="Missingness(%)",x="Sample")+
    ggtitle("All targets.")
  ggsave(paste0(outdir,"/Missingness_boxplot_samplewise.pdf"),width=30,units = "cm")
  
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    filter(orderMean<=nplr_mod.psi) %>%
    ggplot(aes(x=Sample,y=missing_percent))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
    labs(y="Missingness(%)",x="Sample")+
    ggtitle("Putative paralogs removed.")
  ggsave(paste0(outdir,"/Missingness_boxplot_samplewise_paralogsRemoved.pdf"),width=30,units = "cm")
}

plot_paralogy_samplewise<-function(){
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
    labs(y="Paralog Rate (%)",x="Sample")+
    ggtitle("All targets.")
  ggsave(paste0(outdir,"/Paralogy_boxplot_samplewise.pdf"),width=30,units = "cm")
  
  out_meanSort %>% 
    mutate(Sample=fct_reorder(Sample,nplr_mod.pred_resid,mean)) %>% 
    filter(orderMean<=nplr_mod.psi) %>%
    ggplot(aes(x=Sample,y=paralog_percent_ignoreMissing))+
    theme_bw()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
    labs(y="Paralog Rate (%)",x="Sample")+
    ggtitle("Putative paralogs removed.")
  ggsave(paste0(outdir,"/Paralogy_boxplot_samplewise_paralogsRemoved.pdf"),width=30,units = "cm")
  
}
