#R --vanilla --quiet < Morimoto_TDC.2016.Cheng_Vogel_20171125_args.R --args Cheng_Vogel_2016_rep_mean_all_Cheng_Vogel_rep_PC07_.RData             out.setting=F clustering=T ScaleBoot=F  lag.neg.ign=T  genename_protein_description.csv             > Selevsek_2015.all_Qian.log 2>&1 &

# This source code is built under R ver.3.1.2

c_args <- commandArgs(trailingOnly=T)
param.setting     <- c_args[2]
eval(parse(text = param.setting))
#out.setting <- T
if(out.setting==T){
  source(file = './Cheng_Vogel/PG/args/out_setting_Cheng_Vogel.R')
}else{
  c_args <- commandArgs(trailingOnly=T)
  makedataData     <- c_args[1]
  c_clustering     <- c_args[3]
  c_ScaleBoot      <- c_args[4]
  c_lag.neg.ign    <- c_args[5]
  id.convertData   <- c_args[6]
  
  outputDirectry.prefix <- './'
  dataDirectry <- './'
  funcDirectry <- './'
  pkgsDirectry  <- './'
  }

outprefix      <- sub(".RData", "", makedataData)

print(makedataData)
print(c_clustering)
print(c_ScaleBoot)
print(c_lag.neg.ign)
print(id.convertData)


eval(parse(text = c_clustering))
eval(parse(text = c_ScaleBoot))
eval(parse(text = c_lag.neg.ign))

filter.delay <- c(0,1,2)

#```{r setup for HCA and def bootstrap filter }

z_norm            <- TRUE  # Z-Normalize by genes and datanames before E score calc. ?
interTemporal     <- TRUE  # convert data for E score calc. to inter-temporal  ?
rho.interTemporal <- FALSE # convert data for corr.Rho calc. to inter-temporal ?

if(clustering==F){score4HCA='noHCA'}else{
  score4HCA='Qian,2001'
}
methodHclust <- 'ward.D2'
  #'ward.D2'#'correlation'
methodDist   <- 'euclidean'
  #'euclidean'#'average'

alpha.au  =  0.95
alpha.bp  =  0.8
nboot     =  10000
scaling_r <- c(seq(.5,1.5,by=.1))

#```

#```{r prefix for output files}

histogram.pdf.prefix        <- 'histogram'

HCA_by_delay.pdf.prefix     <- 'HCA_delay'
BP_of_edges.csv.prefix      <- 'BP_of_edges_delay'
heatmap_by_delay.pdf.prefix <- 'heatmap_delay' 

line_plot.pdf.prefix        <- 'lines_'
scat_plot.pdf.prefix        <- 'scatt_'

uniProtAccess.csv.prefix       <- 'high_BP_uniProtAccess'

perm_null_rho_histo.pdf.prefix <- 'permNullDistHisto'
permuted_Rho.csv.prefix        <- 'permuted_Rho'

RData_PermRhoAnalysis_save.image.prefix <- 'permuted_Rho'


 # permuted Rho

itt_of_permute <- 10000

#```

#```{r functions sorce code file}

EscoreCalc <- 'func_for_calcEscore_20160123.R' # modified: 2016/01/23
pvclustMod <- 'my.pvclust_20151121.R'
permtRho   <- 'func_for_permutest_of_Rho_20160819.R' 
dataManu   <- 'func_data_mnp_20151125.R'
#```


#```{r Load libraries}

## LIBRARIES
source(file=sprintf('%s/%s',funcDirectry,'require_packages_25NOV17.R'))
detach("package:biomaRt", unload=TRUE)
#```



#Functions for data manupilation
#```{r}
source(file = sprintf('%s%s',funcDirectry,dataManu)) ## added 2015/11/24 SM
## source code was moved to EscoreCalc
#```

#Functions for HCA
#```{r}
source(file = sprintf('%s%s',funcDirectry,EscoreCalc)) ## added 2015/11/24 SM
## source code was moved to EscoreCalc
#```

#Functions for pvclust
#```{r}
source(file = sprintf('%s%s',funcDirectry,pvclustMod)) ## added 2015/11/24 SM
## source code was moved to pvclustMod
# includes 
#```

#Functions for Correlation test
#```{r}
source(file = sprintf('%s%s',funcDirectry,permtRho)) ## added 2015/11/24 SM
## source code was moved to permRho
#```


#```{r Load R_data}

## LOAD DATA


getwd()

load(
  paste(
    dataDirectry,
    makedataData,
    sep='')
  )

#```

#```{r make Z-normalised and time differencial data}

#
# The 'data_ana' is loaded from dataDirectry/makedataData.
# If 'z_norm==TRUE', then the 'data_ana_norm' data is output.
# The 'data_ana_norm is having the value before z-normalisation
# (the variable name is 'val').

# The 'dataGQ' is the data for detect the delay using Qian's method
# in the next chank.
# If 'interTemporal==TRUE' then the value are made from
# 'w.dataGQ <- makeDifData(data_ana_norm %>%...', and this function
# attaches 'd_' before each column names.
#

if(z_norm == TRUE){  # added: 2016/01/23
  data_ana_norm <- data_ana %>%
    rename(val2=val) %>%
    group_by(id,dtname) %>%  
    mutate(val=scale(val2, center = TRUE, scale = TRUE)) %>%
    dplyr::select(-val2)%>%
    ungroup()%>%
    data.frame()
  attributes(data_ana_norm$val) <- NULL
    
  if(interTemporal==TRUE){
    w.dataGQ <- makeDifData(data_ana_norm %>% spread(key=var,value=val),cols) %>%
      dplyr::select(id,dtname,starts_with('d_'))
    timePoint <- timePoint-1
    }
  }

dataGQ <- w.dataGQ[duplicated(w.dataGQ$id)|duplicated(w.dataGQ$id,fromLast=TRUE),]

#```


#```{r classification by delayed time}

w.GQSfled_0 <-  ddply(dataGQ,.(id),gq_method2) %>%
  dplyr::rename(
    delayE0=X1,
    scoreE0=X2,
    delayD0=X1.1,
    scoreD0=X2.1
    )

#```

#```{r histogram }

ftable.delay <- data.frame(
  ftable(w.GQSfled_0$delayE0)
  ) %>%
  mutate(
    delay=as.numeric(as.character(Var1)),
    total=sum(Freq),
    prop=paste(round(Freq/total*100,2),'%',sep=''),
    prop_num = Freq/total*100
  )
print(max(ftable.delay$Freq))
by_y.ax.break <- ifelse(
  round(max(ftable.delay$Freq)/10,-1)==0,
  2,#1,
  round(max(ftable.delay$Freq)/10,0)
  )
y.ax.break <- seq(
  0,
  max(ftable.delay$Freq),
  by_y.ax.break
  )

if(lag.neg.ign==T){  
  # if lag.neg.ign==T then ignor genes whose lag is less than 0
  # and hide percentages (not plot 'plot.type2' layer).
  histo_data <- ggplot(data = ftable.delay %>%
                         filter(delay>=0),
                       aes(x=delay,y=as.numeric(Freq),group=prop))
}else{
  histo_data <- ggplot(data = ftable.delay,
                       aes(x=delay,y=as.numeric(Freq),group=prop))
}

plot.type <- geom_bar(stat="identity")
plot.type2 <- annotate(
  "text", 
  label=ftable.delay$prop, 
  x=ftable.delay$delay, 
  y=ftable.delay$Freq+y.ax.break[2]/2, 
  fontface="italic") 
plot.labs <- labs(x='Timelag',y='Number of genes')
plot.y.axe <- scale_y_continuous(breaks=y.ax.break,limits=c(0,max(ftable.delay$Freq)))
plot.x.axe <- scale_x_continuous(breaks=ftable.delay$delay)
plot.axe.text <- theme(
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=25))#15)) 
pdf(      
  sprintf(
    '%s%s_%s_%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    histogram.pdf.prefix
    ) ,
  width=100,
  paper='USr'
  )
if(lag.neg.ign==T){
  plot( histo_data + plot.type + plot.labs + plot.y.axe + plot.x.axe + plot.axe.text)
}else{
  plot( histo_data + plot.type + plot.labs + plot.y.axe + plot.x.axe + plot.type2 + plot.axe.text)
}
dev.off()

#```
#```{r filtering by time lag}
filter.delay.df <- data.frame(delayE0=filter.delay)

w.GQSfled_0 <- w.GQSfled_0 %>%
  inner_join(filter.delay.df)
GQSfled_0 <- data.frame(ftable(w.GQSfled_0$delayE0)) %>%
  filter(Freq>2)%>%
  mutate(delayE0 = as.numeric(as.character(Var1)))%>%
  dplyr::select(delayE0)%>%
  inner_join(w.GQSfled_0)
  

mst <- dataGQ %>%
  filter(id %in% GQSfled_0$id)

#```

#Hieralchical Clustering Analysis 
#```{r HCA}

#
# The purpose of this large chank is to make the 'picked_genes' data.
# The 'picked_genes' is the filter for selecting genes 

# If 'ScaleBoot==TRUE' then the cluster-wise selection is done
# using Approximately Unbiased and Bootstrap Probability as
# selection criteria.
# Else, all genes are included in the 'picked_genes'.
#

if(ScaleBoot==TRUE){
  sb.switch <- '.on_'
}else{
  sb.switch <- '.off_'
}

if(score4HCA != 'noHCA'){
  if(score4HCA=='Qian,2001'){
    GQScoreMat_0 <-  ddply(mst,.(id),gq_scoreMatOut)#,.progress='text')
    }else
      {
        if(score4HCA %in% c('trc_prt_1dim'))
          GQScoreMat_0 <-  mst %>%
            filter(dtname =='trc') %>%
            inner_join( mst %>%
                          filter(dtname=='prt') %>%
                          setNames(c(names(.)[1], paste0('prt_',names(.)[-1])))
                        ,by='id')%>%
            dplyr::select(-prt_dtname,-dtname)
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
        '%s%s_%s_%s_delay%s_%s_nboot%s_sb%s_output.pdf',
        outputDirectry.prefix,
        outprefix,
        bioproc,
        HCA_by_delay.pdf.prefix,
        delay,
        score4HCA,
        nboot,
        sb.switch
        ),
      width=100,
      paper='USr'
      )
  pv <- silent.pvclust(
    data=t(annMstScoreMat),
    methodHclust,methodDist,nboot
    )
  
  #++++++ scale boot ++ added in 2015/11/27
  pv.sb <- sbfit(pv)
  sink(file = sprintf(fmt =
                      sprintf(
                        '%s%s_%s_%s_delay%s_%s_nboot%s_sb%s_output.txt',
                        outputDirectry.prefix,
                        outprefix,
                        bioproc,
                        HCA_by_delay.pdf.prefix,
                        delay,
                        score4HCA,
                        nboot,
                        sb.switch
                        )
                    )
     )
  print(summary(pv.sb,k=c(1:3)))
  sink()
  
  if(ScaleBoot==TRUE){
    pv  <- sbpvclust(pv, pv.sb)
    }
  #++++++ end (added in 2015/11/27)

  plot(pv,cex.pv=0.1,lwd=0.5,cex=0.1)
  for(edges in 1:nrow(pv$edges)){
     msplot(pv,edges)
    }
  pvcl_se <- seplot(pv, identify=TRUE)
  plot(pv,cex=0.5)
  dev.off()

    pv.picked <- my.pvrect(
      pv,
      alpha.au = alpha.au ,
      alpha.bp = alpha.bp ,
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
      '%s%s_%s_%s_delay%s_%s_nboot%s_sb%s.csv',
      outputDirectry.prefix,
      outprefix,
      bioproc,
      BP_of_edges.csv.prefix,
      delay,
      score4HCA,
      nboot,
      sb.switch
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
      gq_scoreMatOut_withName) %>% 
    llply(my.heatmap.3)
  dev.off()  
  }
  picked_genes <- picked_genes %>%
    dplyr::rename(id=pv.picked.clusters..i...j., edges=pv.picked.edges..i..) %>%
    filter(delay %in% filter.delay)
  
  }else{ # end:if(score4HCA != 'noHCA')
    picked_genes <- data.frame()
    pv.picked <- GQSfled_0 %>%
      dplyr::select(id,delayE0) %>%
      dplyr::rename(delay=delayE0) %>%
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
      dplyr::rename(
        id=pv.picked.id..i.., 
        edges=pv.picked.edges..i..,
        delay=pv.picked.delay..i..) %>%
      filter(delay %in% filter.delay)
    }

dataPath_1     <- 'HighAUcluster_20171208.csv' # same definition as 'biomaRt_query_for_...'

if(score4HCA != 'noHCA'){
  HighAUcluster <- picked_genes
  write.csv(
    HighAUcluster,
    sprintf(
      '%s/%s', 
      dataDirectry,
      dataPath_1
      )
    )
}


#```

#```{r selecting genes using picked.genes data}

picked_genes <- read.csv(
  sprintf(
    '%s/%s', 
    dataDirectry,
    dataPath_1
    )
  )

source(
  sprintf(
    '%s/%s',
    funcDirectry,
    '_Cheng_Vogel_biomaRt_query_for_20171125_output.R'
    )
  )

inputName_to_outputName <- picked_genes %>%
  inner_join(
    read.csv(
      sprintf('%s%s',
              dataDirectry,
              id.convertData
              )
      )
    ) %>%
  dplyr::rename(accession=external_gene_name.x)

edge_to_inputName_to_outputName = picked_genes %>%
  full_join(inputName_to_outputName)

sum(is.na(edge_to_inputName_to_outputName))
sum(duplicated(edge_to_inputName_to_outputName$id))

pander(edge_to_inputName_to_outputName)
write.csv(
  edge_to_inputName_to_outputName,
  sprintf(
    '%s%s_%s_%s_%s_nboot%s_sb%s.csv',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    uniProtAccess.csv.prefix,
    score4HCA,
    nboot,
    sb.switch
    )
  )

#```


#```{r Spearmans Rho }

# If inputName has duplicated outputName, 
# then choose the first one.

df.genes <- ddply( 
  edge_to_inputName_to_outputName[
    !(duplicated(edge_to_inputName_to_outputName$id)),],
  .(delay,edges),
  paste
  ) %>%
  dplyr::select(-V2,-V3)

df.genes.notDelay <- ddply( 
  edge_to_inputName_to_outputName[
    !(duplicated(edge_to_inputName_to_outputName$id)),],
  .(edges),
  paste
  ) %>%
  dplyr::select(-V2,-V3)

#```{r making annMstLong data}
#
# The 'edge.dataGQ' is the data for plots and Rho calculation.
# If 'z_norm==TRUE', then the 'dataGQ' is "makeDifData(data_ana_norm %>% ...".

# And, if 'rho.interTemporal == TRUE', then the 'edge.dataGQ' is
# made from interTemporal data of z-normalised data.

# If 'rho.interTemporal != TRUE', then the 'edge.dataGQ' is made from
# the 'data_ana', not-normalised data. Then, the 'annMstLong' is
# made from this 'edge.dataGQ' with normalised values
# ($norm.trc, $norm.prt) having not-normalised values($trc, $prt).

# The normalised data in 'data_ana_norm' is normalised across the all
# timepoints, BUT in 'annMstLong' across the selected timepoints 
# (in 'cols').
#

if(rho.interTemporal==TRUE){
  edge.dataGQ <-
    inner_join(
      dataGQ,
      edge_to_inputName_to_outputName,
      by='id'
    )
  annMstLong <- edge.dataGQ %>%
    gather(var,val,starts_with('d_')) %>%
    spread(dtname,val) %>%
    mutate(
      var.prt = match(
        var,
        names(
          edge.dataGQ %>%
            dplyr::select(
              starts_with('d_')
            )
        )
      ), 
      var.trc = var.prt + delay
    )
  }else{
    
    edge.dataGQ <-
      inner_join(
        data_ana %>%
          filter(var %in% cols) %>%
          mutate(var=factor(var,cols)) %>% 
          spread(key=var,value=val),
        edge_to_inputName_to_outputName,
        by='id'
        )
    
    annMstLong <- edge.dataGQ %>%
      gather(var,val,one_of(cols)) %>%
      spread(dtname,val) %>%
      group_by(id) %>%
      mutate(
        var.prt = match(
          var,
          names(
            edge.dataGQ %>%
              dplyr::select(
                one_of(cols)
                )
            )
          ), 
        var.trc = var.prt + delay,
        prt.norm=scale(prt,center=TRUE,scale=TRUE),
        trc.norm=scale(trc,center=TRUE,scale=TRUE),
        var     =factor(var, cols)
        )%>%
      ungroup()
    attributes(annMstLong$prt.norm) <-NULL
    attributes(annMstLong$trc.norm) <-NULL
    }


delayAnnMstLong <- inner_join(
  annMstLong %>% dplyr::rename(var.new=var.prt) %>% 
    dplyr::select(id,prt,var.new,prt.norm) ,
  annMstLong %>% dplyr::rename(var.new=var.trc) %>%
    dplyr::select(id,trc,var.new,trc.norm,var,edges,delay,accession) ,
  by = c("id", "var.new")
  ) %>%
  dplyr::rename(var.ori=var)

write.csv(
  x = data.frame(annMstLong),
  file = sprintf(
    fmt = '%s%s_data.csv',
    outputDirectry.prefix,outprefix)
    )
#```

#```{r plotting line plot}

# make pdf file of time series line plots

line.plot <- function(data){
  id <- data$id[1]
  accession <- data$accession[1]
  delay <- data$delay[1]
  edge  <- data$edges[1]
  p <- ggplot(data %>%
                mutate(timepoint=paste('timepoint #',var.ori)),
              aes(x=var.ori,y=val,group=var)
              )
  plot.type <- geom_line(aes(colour=var),size=1.5)
  plot.ylim <-coord_cartesian(ylim=c(-2.5,2.5))
  plot.theme <- theme(
    legend.background=element_blank(),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(colour="black", fill="white"),
    axis.text.x = element_text(size=30),
    axis.text.y = element_text(size=25)
    )
  plot.bg    <- theme_bw()
  plot.title <- labs(
    list(
      title = paste(
        'Gene name_',id,'(',accession,');','Time-lag ; ',delay,'[Cluster;',edge,']',
        sep=''
        )
      )
    )
  plot.ylab <- ylab('Normalised')
  plot( p+ plot.type+ plot.theme + plot.title + plot.ylim + plot.ylab)
  }

pdf(
  sprintf(
    '%s%s_fixed_y_%s_%s_shifted_%s_nboot%s_sb%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    line_plot.pdf.prefix,
    score4HCA,
    nboot,
    sb.switch
    ),
  width=16,
  height=10
  )
  line.plot.list <- 
    dlply(
      delayAnnMstLong %>% 
        gather(var,val,prt.norm,trc.norm) ,
      .(delay,edges,id),
      line.plot
      ) 
dev.off()

pdf(
  sprintf(
    '%s%s_fixed_y_%s_%s%s_nboot%s_sb%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    line_plot.pdf.prefix,
    score4HCA,
    nboot,
    sb.switch
  ),
  width=16,
  height=10
)
line.plot.list <- 
  dlply(
    annMstLong %>%
      rename(var.ori=var)%>%
      gather(var,val,prt.norm,trc.norm) ,
    .(delay,edges,id),
    line.plot
  ) 
dev.off()


#```

#```{r plotting scatter plot}

scat.plot <- function(data){
  shape_basket <- c(8,15,16,17,18,19)
#  shape_basket <- c(9:14) 
  id <- data$id[1]
  accession <- data$accession[1]
  delay <- data$delay[1]
  edge  <- data$edges[1]
  RHO   <- cor.test_by.id(data,'prt','trc')$rho.est
  p <- ggplot(data %>%
                mutate(
                  timepoint=paste('timepoint #',var.new),
                  cluster=paste('Cluster ;',edges),
                  delay=paste('Time-lag ; ',delay,sep='')) ,
              
              aes(x=trc,y=prt,group=timepoint))
  
  if(length(unique(data$id))<=6){
    plot.type <- geom_point(aes(colour=timepoint,shape=id),size=20,alpha=0.6,stroke=5.0)
  }else{
    plot.type <- geom_point(aes(colour=timepoint,shape=id),size=8,alpha=0.8)
  }

  if(length(unique(data$id))<=6){
    plot.shape <- scale_shape_manual(
      values=shape_basket[c(1:length(unique(data$id)))]
      )
  }else{
    plot.shape <- scale_shape_manual(
      values=rep(16,length(unique(data$id))))
  }
  
  legends.guide <-guides(
    guide_legend(keywidth = 0.3, keyheight = 0.3,nrow = 2, byrow = TRUE)
    )
  plot.ylim <-coord_cartesian(ylim=c(-2.5,2.5))
  plot.facet <- facet_wrap(~edges+delay)
  plot.title <- labs(
    list(
      title = paste(
        'Rank Correlation (Rho;',
        round(RHO,digits = 4),
        ')'
        )
      )
    )
  plot.bg <- theme_gray()
  plot.axe.text <- theme(
    axis.text.x = element_text(size=30),
    axis.text.y = element_text(size=30)) 
  plot( p+ plot.type+ plot.shape + plot.ylim + plot.facet+ plot.title + plot.bg +legends.guide + plot.axe.text )
}

pdf(
  sprintf(
    '%s%s_%s_%s_%s_nboot%s_sb%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    scat_plot.pdf.prefix,
    score4HCA,
    nboot,
    sb.switch
    )
  ,width=16,
  height=10
  )
  scat.plot.list <- dlply(
    data_ana_norm %>%        # normalised in 'if(z_norm==TRUE){...}'
      filter(dtname=='prt') %>% rename(prt=val) %>%
      inner_join(
        data_ana_norm %>% 
          filter(dtname=='trc') %>% rename(trc=val),
        by=c('id','var')
      ) %>% 
      mutate(
        delay='disabled',
        edges='All Genes',
        var.new=factor(var)
      ) %>%
      
      # --- added for _Cheng_Vogel_analysis --- #
      
      inner_join(
        inputName_to_outputName  %>%
          dplyr::select(id, accession),
        by='id'
      )  %>% 
      dplyr::select(-id) %>%
      dplyr::rename(id = accession)
    
     # ---------------------------------------- #
    ,
    .(delay), 
    scat.plot) # ALL GENES (normalised in 'if(z_norm==TRUE){...}')
  
  scat.plot.list <- dlply(
    annMstLong %>%        # normalised after 'picked.genes' selection.
      rename(raw.prt=prt,raw.trc=trc)%>%
      rename(prt=prt.norm,trc=trc.norm)%>%
      mutate(
        delay='disabled(genes selected as analyte)',
        var.new=factor(var))%>%
      
      # --- added for _Cheng_Vogel_analysis --- #
      dplyr::select(-id) %>%
      dplyr::rename(id = accession),
    # ---------------------------------------- #
    .(delay), 
    scat.plot
    ) #
  
  scat.plot.list <- dlply(
    annMstLong %>%
      rename(raw.prt=prt,raw.trc=trc)%>%
      rename(prt=prt.norm,trc=trc.norm)%>%
      mutate(delay='disabled',var.new=factor(var))%>%
      
      # --- added for _Cheng_Vogel_analysis --- #
      dplyr::select(-id) %>%
      dplyr::rename(id = accession),
    # ---------------------------------------- #
    .(edges),
    scat.plot
    ) 
  
  scat.plot.list <- dlply(
    delayAnnMstLong %>%
      dplyr::rename(raw.prt=prt,raw.trc=trc)%>%
      dplyr::rename(prt=prt.norm,trc=trc.norm)%>%
      mutate(var.new=factor(var.ori))%>%
      
      # --- added for _Cheng_Vogel_analysis --- #
      dplyr::select(-id) %>%
      dplyr::rename(id = accession),
    # ---------------------------------------- #
      .(delay,edges), scat.plot)
  
dev.off()

rho.df <- ddply(delayAnnMstLong,.(edges,delay),
                cor.test_by.id,
                'trc.norm','prt.norm',method='s'
                )%>%
  inner_join(
    df.genes %>% dplyr::select(V1,V4,delay,edges) ,
    by=c('delay','edges')
    )
rho.df.notDelay <- ddply(annMstLong %>% 
                           mutate(delay='disable',var.new=factor(var)),
                         .(edges),
                cor.test_by.id,
                'trc.norm','prt.norm',method='s'
                )%>%
  inner_join(
    df.genes.notDelay %>% dplyr::select(V1,V4,edges) ,
    by=c('edges')
    )
#```

#Empirical test of Rho /n

#NULL distribution of each cluster is created /n
#by the frequency of Rho from `r itt_of_permute` permutations)

#```{r empirical test for Rho using permutation test  }

# make data for permutest

empi.rho_id_Shufle <- function(data){
  for(i in 1:itt_of_permute){
    if(i!=1){
      mstSfl <- ddply(
        data,.(id,dtname),
        varSfl,sfl=TRUE,startCol=3,endCol=length(mst))#,.progress='text')
      } else {
        mstSfl <- data
    }
    if (rho.interTemporal==TRUE){
      mstSflLong <- mstSfl %>%
        gather(var,val,starts_with('d_')) %>%
        spread(dtname,val) %>%
        mutate(
          var.prt = match( var, 
                           names(
                             edge.dataGQ %>%
                               dplyr::select(
                                 starts_with('d_')
                               )
                           )
          ),
          var.trc = var.prt + delay
        )
    }else{
      mstSflLong <- mstSfl %>%
        gather(var,val,one_of(cols)) %>%
        spread(dtname,val) %>%
        mutate(
          var.prt = match( var, 
                           names(
                             edge.dataGQ %>%
                               dplyr::select(
                                 one_of(cols)
                               )
                           )
          ),
          var.trc = var.prt + delay
        )
    }

    
    
    delayMstSflLong <- inner_join( 
      # inner_join deletes protrusion time point
      mstSflLong %>% rename(var.new=var.prt) %>%
        dplyr::select(id,prt,var.new) ,
      mstSflLong %>% rename(var.new=var.trc) %>%
        dplyr::select(id,trc,var.new) ,
      
      by=c('id','var.new')
      )
  
  GQSfled <-  ddply(
    delayMstSflLong,.(),cor.test_by.id,'prt','trc'
    )
  
    #,.progress='text')
  GQSfled$itt <- i
  
  if(i==1){
    GQSfled_i <- GQSfled
  }else{
    GQSfled_i <- rbind(GQSfled,GQSfled_i)
    }
  #print(i)
  }
  #for
return(GQSfled_i)
}


output.empi.rho <- ddply(edge.dataGQ %>%
                           gather(var,val,one_of(cols))%>%
                           group_by(id,dtname)%>%
                           mutate(
                             norm=scale(val,scale=TRUE,center=TRUE) 
                             )%>%
                           ungroup()%>%
                           data.frame() %>%
                           dplyr::select(-val)%>%
                           spread(key = var,value = norm) %>%
                           dplyr::select(
                             id,
                             dtname,
                             one_of(cols),
                             edges,
                             delay,
                             accession
                             ),
                         .(delay,edges),
                         empi.rho_id_Shufle
                         ) %>%
  inner_join(df.genes,by=c('delay','edges')) %>%
  ddply(.(delay,edges),my.rank,'rho.est') %>%
  mutate(
    p_value = 1-(rank / itt_of_permute)
    ) 

panderOptions('table.split.table', 1000)
panderOptions('table.continues', '')
pander(
  output.empi.rho %>%
    filter(itt==1 ) %>%
    dplyr::select(rho.est,p_value,V1,V4,delay,edges)
  )

output.empi.rho.q <- output.empi.rho %>%
  filter(itt==1 ) %>%
  mutate(q_value=p.adjust(p_value,method = 'BH'))

write.csv(output.empi.rho.q,
          file=sprintf(
            '%s%s_%s_%s_%s_nboot%s_sb%s.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            permuted_Rho.csv.prefix,
            score4HCA,
            nboot,
            sb.switch
            ) 
          )

#```



#```{r Histogram}
## Histogram

head(output.empi.rho)
qtl.empi.rho <- output.empi.rho %>%
  ddply(.(delay,edges),my.qtl,'rho.est') %>%
  rename(
    Q.025=V1,
    Q.050=V2,
    Q.950=V3,
    Q.975=V4
  ) %>%
  gather(var,val,Q.950)#starts_with('Q'))

df.perm_null_rho_histo <- inner_join(qtl.empi.rho,output.empi.rho)

perm_null_rho_histo <- function( data){
  
  delay <- unique(data$delay)
  edges <- unique(data$edges)
  print(c(delay,edges))
  
  p     <- ggplot(data ,aes(x=rho.est))
  p     <- p + geom_histogram(position='identity',binwidth = 0.01)
  
  gVline   <- geom_vline(
    aes(
      xintercept = rho.est,
      colour     = 'red',
      linetype='realised'),
    data = data %>%
      filter(itt==1),
    show_guide=T
  ) 
  
  qtlVline <- geom_vline( 
    aes(
      xintercept = val,
      linetype=var
    ),
    data=data[,c('var','val')] ,
    show_guide=T
  ) 
  #  facet    <- facet_wrap(~id,scales='free')
  
  legend.title <- labs(
    title=paste('null dist. of Rhos from ',
                itt_of_permute,
                ' permutations; delay=',delay,',Edge=',edges,sep=''),
    linetype='vert.lines'
  )
  legend.theme <- theme(
    legend.text = element_text(size = 20, colour = "black", angle = 45),
    strip.text.x = element_text(size =20, colour = "black", angle = 0),
    legend.background=element_blank(),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(colour="black", fill="white"))
  
  gg.result <- p + gVline + qtlVline + legend.theme + legend.title
  plot(gg.result)
  
}
pdf(
  sprintf(
    '%s%s_%s_%s_%s_nboot%s_sb%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    perm_null_rho_histo.pdf.prefix,
    score4HCA,
    nboot,
    sb.switch
  )
)
perm.rho.histo <- dlply(
  df.perm_null_rho_histo,
  .(delay,edges),
  perm_null_rho_histo
)
dev.off() ## pdf(perm_null_rho_histo)

#```
read.csv(
  sprintf(
    '%s%s_%s_%s_%s_nboot%s_sb%s.csv',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    uniProtAccess.csv.prefix,
    score4HCA,
    nboot,
    sb.switch
    ) 
  ) %>%
  inner_join(
    read.csv(
      file=sprintf(
        '%s%s_%s_%s_%s_nboot%s_sb%s.csv',
        outputDirectry.prefix,
        outprefix,
        bioproc,
        permuted_Rho.csv.prefix,
        score4HCA,
        nboot,
        sb.switch
        )
      ) %>%
      dplyr::select(delay, edges, rho.est, q_value) %>%
      mutate(
        rho = round(
          rho.est, 2),
        q_value = ifelse(
          q_value==0, 
          format(round(0.0001,4),scientific=FALSE), 
          format(round(q_value,4),scientific=FALSE)
          )
        )
    ) %>%
  mutate( 
    analysis_res = paste( rho, '\r', '[', q_value, ']', sep=''),
    description = sub(
      ':',
      ':\r',
      description
      )
    ) %>%
  mutate(
    description = sub(
      "\\[Source:",
      "\r\\[Source:",
      description
      )
    ) %>%
  dplyr::select( id, edges, delay, accession, description, q_value, rho, analysis_res) %>%
  write.csv(
    sprintf(
      '%s_Analyses_results_%s_%s_%s_%s_nboot%s_sb%s.csv',
      outputDirectry.prefix,
      outprefix,
      bioproc,
      uniProtAccess.csv.prefix,
      score4HCA,
      nboot,
      sb.switch
      )
  )
      
save.image(
  sprintf(
    '%s_AnalysisResults_%s.RData',
    outputDirectry.prefix,
    outprefix
    )
  )
