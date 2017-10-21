
#=================================================================================
#   filt_high.AU.BP

# input *delay*,  

#=================================================================================

filt_high.AU.BP <- function(score4HCA){
  if(score4HCA != 'noHCA'){
    if(score4HCA=='Qian,2001'){
      GQScoreMat_0 <-  ddply(mst,.(id),gq_scoreMatOut)#,.progress='text')
    }else
    {
      if(score4HCA %in% c('trc_prt_1dim', 'noHCA'))
        GQScoreMat_0 <-  mst %>%
          filter(dtname =='trc') %>%
          inner_join( mst %>%
                        filter(dtname=='prt') %>%
                        setNames(c(names(.)[1], paste0('prt_',names(.)[-1])))
                      ,by='id')%>%
          dplyr::select(-prt_dtname)
    }
    all_GQSfled_0_delayed <- GQSfled_0 %>%   
      filter(delayE0>=0)
    group_by_delay <- unique(all_GQSfled_0_delayed$delayE0)
    picked_genes <- data.frame()
    
    for (delay in group_by_delay){
      GQSfled_0_delayed <- 
        all_GQSfled_0_delayed[order(all_GQSfled_0_delayed$scoreE0),] %>%
        filter(delayE0==delay)
      mstScoreMat <- GQScoreMat_0 %>%
        filter(id %in% GQSfled_0_delayed$id)  
      dimnames(mstScoreMat)[[1]] <- mstScoreMat$id
      annMstScoreMat <- mstScoreMat[,-1]
      
      HCA_Ward <- hclust(
        dist(annMstScoreMat),
        method=methodHclust
      )
      plot(HCA_Ward,cex = 0.5)
      
      pdf(
        sprintf(
          '%s%s_%s_%s_delay%s_%s_nboot%s_output.pdf',
          outputDirectry.prefix,
          outprefix,
          bioproc,
          HCA_by_delay.pdf.prefix,
          delay,
          score4HCA,
          nboot
        ),
        width=100,
        paper='USr'
      )
      pv <- pvclust(
        t(annMstScoreMat),
        method.hclust=methodHclust,
        method.dist=methodDist,
        nboot=nboot,
        r=multScale,
        quiet=TRUE
        ) 
      plot(pv,cex.pv=0.1,lwd=0.5,cex=0.1)
      for(edges in 1:nrow(pv$edges)){
        msplot(pv,edges)
      }
      
      # p-value vs standard error plot
      pvcl_se <- seplot(pv, identify=TRUE)
      plot(pv,cex=0.5)
      dev.off()
      
      pv.picked <- my.pvrect(
        pv,
        alpha.au = 0.9 ,
        alpha.bp = 0.8 ,
        pv='bp',
        type='gt'
      )
      if(!is.null(pv.picked$clusters[[1]])){
        for( i in 1:length(pv.picked$edges)){
          for( j in 1:length(pv.picked$clusters[[i]])){
            d.f_ij <-  data.frame(pv.picked$clusters[[i]][j],
                                  pv.picked$edges[[i]],
                                  delay
            )
            picked_genes <- rbind(picked_genes,d.f_ij)
          }
        }
      }
      
      pv_df <- data.frame(pv$edges)
      write.csv(
        pv_df, 
        file=sprintf(
          '%s%s_%s_%s_delay%s_%s_nboot%s_output.csv',
          outputDirectry.prefix,
          outprefix,
          bioproc,
          BP_of_edges.csv.prefix,
          delay,
          score4HCA,
          nboot
        )
        , col.names=T
        , quote=F
        , row.names=T
      )
      
      
      pdf(sprintf('%s%s_%s_%s_%s_output.pdf',
                  outputDirectry.prefix,
                  outprefix,
                  bioproc,
                  heatmap_by_delay.pdf.prefix,
                  delay
      ),
      width=8,
      height=8
      )
      GQScoreHeatmap_0 <-
        dlply(
          mst %>% filter(id %in% dimnames(annMstScoreMat)[[1]]), 
          .(id),
          gq_scoreMatOut) %>% #,.progress='text')
        llply(matrix,nrow=timePoint+1,byrow=TRUE) %>%
        llply(heatmap.3,trace="none", dendrogram="none", Rowv=F, Colv=F, color.FUN="redgreen", 
              cluster.by.row=F, cluster.by.col=F, mapratio=1, mapsize=4, main='' )
      dev.off()  
    }
    picked_genes <- picked_genes %>%
      rename(id=pv.picked.clusters..i...j., edges=pv.picked.edges..i..) %>%
      filter(delay %in% filter.delay)
    
  }else{ # end:if(score4HCA != 'noHCA')
    picked_genes <- data.frame()
    pv.picked <- GQSfled_0 %>%
      dplyr::select(id,delayE0) %>%
      rename(delay=delayE0) %>%
      mutate(edges=bioproc) %>%
      dplyr::select(id,edges,delay)
    for( i in 1:length(pv.picked$edges)){
      d.f_ij <-  data.frame(pv.picked$id[[i]],
                            pv.picked$edges[[i]],
                            pv.picked$delay[[i]]
      )
      picked_genes <- rbind(picked_genes,d.f_ij)
    }
    picked_genes <- picked_genes %>%
      rename(
        id=pv.picked.id..i.., 
        edges=pv.picked.edges..i..,
        delay=pv.picked.delay..i..) %>%
      filter(delay %in% filter.delay)
  }
}


