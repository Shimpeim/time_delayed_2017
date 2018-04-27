#R --vanilla --quiet < simu.permtest_test_20160817_args.R --args  50 10000 YDL110C_YDR122W Selevsek_2015_all_Selevsek_20160706.RData batch_1  > Simulation_permtest_20160817_1.log 2>&1 &

# save.image(
#     sprintf(
#        '%s/%s',
#        outputDirectry.prefix,
#        'simu.permtest_test_20160803_.RData'
#        )
#     )

c_args <- commandArgs(trailingOnly=T)
c_sampleSize       <- c_args[1]
c_itt_of_permute   <- c_args[2]
c_id.select        <- c_args[3] # FMT : geneName1_geneName2_ ...
c_makedataData     <- c_args[4]
c_outdat_surfix    <- c_args[5]

#
#
#c_sampleSize       <- '100'
#c_itt_of_permute   <- '10000'
#c_id.select        <- 'a_b' #_c' #_d_e_f_g_h_i_j_k_l_m_n_o_p_q_r_s_t_u_v_w_x_y_z'
#c_makedataData     <- 'Selevsek_2015_all_Selevsek_20160706.RData'
#
#

sim.cols <- paste('V',c(1:5),sep='')
itt_of_permute <- as.numeric(c_itt_of_permute)
makedataData <- c_makedataData

sampleSize <- as.numeric(c_sampleSize) # N of random sampled data for each rho 

outprefix      <- sprintf('%s_permute_%s_%s','alpha_error_simulation_for_CV2016_permtest_itt_',c_itt_of_permute,c_outdat_surfix)
bioproc        <- sprintf('%s','_MultValNorm_cov_')
outputDirectry.prefix          <- './'
#outputDirectry.prefix         <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test/test/20160817'

dataDirectry <- './'
funcDirectry <- './'
pkgsDirectry  <- './'

#dataDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test/test/20160817/'
#funcDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test/test/20160817/'
#pkgsDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test/'


rho.interTemporal <- FALSE # convert data for corr.Rho calc. to inter-temporal ?

#```{r prefix for output files}

histogram.pdf.prefix        <- 'histogram'

scat_plot.pdf.prefix        <- 'scatt_'
line_plot.pdf.prefix        <- 'line_'

perm_null_rho_histo.pdf.prefix <- 'permNullDistHisto'
permuted_Rho.csv.prefix        <- 'permuted_Rho'

RData_PermRhoAnalysis_save.image.prefix <- 'permuted_Rho'

#```

#```{r functions sorce code file}

permtRho   <- 'func_for_permutest_of_Rho_20160819.R' 
dataManu   <- 'OldFunc_20151125.R'
EscoreCalc <- 'func_for_calcEscore_20160123.R' # modified: 2016/01/23
#```


#```{r Load libraries}

## LIBRARIES

packages <- c(
  'MASS',
  'dplyr', # progress bar
  'plyr',  # progress bar
  'dplyr', # progress bar
  'tidyr',
  'stringr',
  'xlsx',
  'ggplot2',
  'gplots',
  'GMD',
  'reshape2',
  'pander',
  'matrixcalc'
) 

new.packages <-
  packages[!(packages %in% installed.packages())] 
# installed.packages() returns installed packages 

if(length(new.packages) > 0){ 
  install.packages(new.packages, repos='http://cran.us.r-project.org')
}
require('MASS')
require('plyr')  # progress bar
require('dplyr') # progress bar
require('tidyr')
require('stringr')
require('xlsx')
require('ggplot2')
require('gplots')
require('GMD')
require('pvclust')
require('reshape2')
require('pander')
require('matrixcalc')

#```

id.select <- unlist(                   # Genes in the analysing gene-cluster
  str_split(c_id.select, '_', n=Inf)
)

#Functions for HCA
#```{r}
source(file = sprintf('%s%s',funcDirectry,EscoreCalc)) ## added 2015/11/24 SM
## source code was moved to EscoreCalc
#```

#Functions for data manupilation
#```{r}
source(file = sprintf('%s%s',funcDirectry,dataManu)) ## added 2015/11/24 SM
## source code was moved to EscoreCalc
#```

#Functions for Correlation test
#```{r}
source(file = sprintf('%s%s',funcDirectry,permtRho)) ## added 2015/11/24 SM
## source code was moved to permRho
#```

##====================##
## Simulation dataset ##
##====================##

# Sigma matrix from Observed data

load(
  file = sprintf(
    '%s%s',
    dataDirectry,
    makedataData
    )
  )

data_ana_norm <- dataGQ[
    duplicated(dataGQ$id)|duplicated(dataGQ$id,fromLast=TRUE),
    ] %>%
  gather(var,val,-id,-dtname) %>%
  rename(val2=val) %>%
  group_by(id,dtname) %>%  
  mutate(val=scale(val2, center = TRUE, scale = TRUE)) %>%
  dplyr::select(-val2)%>%
  ungroup()%>%
  data.frame() %>%
  spread(key = var,value=val)

norm_trc <- data_ana_norm%>%
  filter(dtname=='trc')
norm_prt <- data_ana_norm%>%
  filter(dtname=='prt')

Sigma_norm_trc  <- cov(norm_trc[,cols])#,method='s')
Sigma_norm_prt  <- cov(norm_prt[,cols])#,method='s')



#

SimSet <- sampleSize
id.select <- unlist(                   # Genes in the analysing gene-cluster
  str_split(c_id.select, '_', n=Inf)
  )

cross.corr.mvr_obscor <- function(output.name,SimSet,dtname,Sigma_mat){
  X_sigma_set  <- data.frame()
  SimSet_update <- SimSet
  count_update  <- 0
  sim.cols <- cols
  Sigma <- Sigma_mat
  

    X_all <- as.data.frame(
      mvrnorm(n = SimSet_update*length(id.select), 
              rep(0, length(sim.cols)),
              Sigma_mat,
              tol = 1e-10,
              empirical=T
      )
    ) %>%
      cbind(
        count    =rep(
          seq( from=count_update+1,to=sum(count_update,SimSet),by=1 ),
          each=length(id.select)
        ),
        gene   = rep(id.select,SimSet),
        dtname = rep(
          c(dtname),
          SimSet*length(id.select))
      ) %>%
      gather(
        key = var,
        value=val,
        -count,-dtname,-gene
      ) %>%
      mutate(
        var.ori = match(var,sim.cols)  , 
        var.new = match(var,sim.cols) + 1
      )
    
    X_sigma <- X_all %>%
      inner_join(
        X_all,
        by=c(
          'var.ori'='var.new',
          'count','dtname','gene'
          )
      ) %>%
      dplyr::select(
        -var.ori,-var.ori.y,-var.new
      ) %>%
      ddply(
        .(count,dtname),
        cor.test_by.id,
        'val.y','val.x',
        method.cor='spearman',
        .progress='text'
      ) %>%
      inner_join(X_all,by=c('count','dtname')) 
    X_sigma_set    <- rbind(X_sigma_set,X_sigma)

  return(X_sigma_set)
}

sim_norm_10000_trc <- cross.corr.mvr_obscor(
  sim_norm_10000,
  SimSet=sampleSize,
  dtname='trc',
  Sigma_mat=Sigma_norm_trc
  ) %>%
  rename(count2 = count)

sim_norm_10000_prt <- cross.corr.mvr_obscor(
  sim_norm_10000,
  SimSet=sampleSize,
  dtname='prt',
  Sigma_mat=Sigma_norm_prt
  )%>%
  rename(count2 = count)

######################## 2016/8/18
########################

write.csv(
  rbind(
    sim_norm_10000_trc,
    sim_norm_10000_prt
    ),
          file = sprintf(
            '%s%s_%s_%s.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            'simu_dat_20160818')
)


##================================##
## Choose 10 representation data  ##
##================================##

sim.dat.sample <- rbind(
  data.frame(  sim_norm_10000_trc,rho.which='observed'),
  data.frame(  sim_norm_10000_prt,rho.which='observed')
  ) %>%
  filter(
    count2 %in% c(1:10)
    )  %>%
  mutate(id = paste('SimNo.',count2,sep='')) %>%
  dplyr::select(-count2) 
  

sim.dat.sample_rho.est <- sim.dat.sample %>%
  filter(dtname=='trc') %>%
  inner_join(sim.dat.sample %>%
               filter(dtname=='prt'),
             by=c(
               'id',
               "var",
               "var.ori",
               "var.new",
               'gene',
               "rho.which"
               )
             ) %>%
  ddply(
    .(id,rho.which),
    cor.test_by.id,
    'val.x',
    'val.y'
    ) %>%
  rename(
    rho= rho.est
    ) %>%
  inner_join(
    sim.dat.sample,
    by=c('rho.which','id')
    )

#```{r make pdf file of time series line plots}
line.plot <- function(data){
  theme_set(theme_classic(base_size = 18, base_family = 'Helvetica'))
  id <- data$id[1]
  rho.which <- data$rho.which
  rho       <- round(data$rho,2)
#  rho.bwtp  <- round(data$rho.bw.tp,2)
  rho.bwtp  <- ''
  p <- ggplot(data %>%
                mutate(timepoint=paste('timepoint ',var)),
              aes(x=var,y=val,group=dtname)
  )
  plot.type <- geom_line(aes(colour=dtname),size=1.5)
  plot.theme <- theme(
    legend.background=element_blank(),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(colour="black", fill="white"))
  plot.bg    <- theme_bw()
  plot.title <- labs(
    list(
      title = paste(
        id,' [Rho(intertemporal)=',rho,'(',rho.bwtp,')]',#,count:',edge,']',
        sep=''
      )
    )
  )
  plot( p+ plot.type+ plot.theme + plot.title)
}

pdf(
  sprintf(
    '%s%s_%s_%s_.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    line_plot.pdf.prefix
  ),
  width=16,
  height=10
)
line.plot.list <- 
  dlply(
    sim.dat.sample_rho.est,
    .(rho.which,id,gene),
    line.plot
  ) 
dev.off()


#```{r plotting scatter plot}

scat.plot_permu.simu <- function(data,x.lab,y.lab){

  shape_basket <- c(8,15,16,17,18,19) # plot shape selecter. (ref. help(shape) {ggplot2})

  id   <- as.character(data$id[1])
  gene <- as.character(data$gene)
  rho.which <- as.character(data$rho.which[1])
  rho       <- round(data$rho[1],2)

  p <- ggplot(data %>%
                mutate(timepoint=paste('timepoint #',var)
                ),
              aes(x=x,y=y,group=timepoint))
  
  if(length(unique(data$id))<=6){
    plot.type <- geom_point(aes(colour=timepoint,shape=gene),size=15,alpha=0.6)
  }else{
    plot.type <- geom_point(aes(colour=timepoint,shape=gene),size=4,alpha=0.6)
  }
  
  if(length(unique(data$gene))<=6){
    plot.shape <- scale_shape_manual(values=shape_basket[c(1:length(unique(data$gene)))])
  }else{
    plot.shape <- scale_shape_manual(values=rep(18,length(unique(data$gene))))
  }
  
  legends.guide <-guides(
    guide_legend(keywidth = 0.3, keyheight = 0.3,nrow = 2, byrow = TRUE)
  )
  plot.title <- labs(
    list(
      title = paste(
        id," [Spearman's r(intertemporal)=",rho,'(',rho.which,')]',#,count:',edge,']',
        sep=''
      ),
      x=x.lab,y=y.lab
    )
  )
  plot.bg <- theme_bw()
  plot( p+ plot.type+ plot.shape + plot.title + plot.bg +legends.guide )
}

pdf(
  sprintf(
    '%s%s_%s_%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    scat_plot.pdf.prefix
  )
  ,width=16,
  height=10
)
scat.plot.list <- dlply(
  sim.dat.sample_rho.est %>% 
    dplyr::select(
      val, id, dtname, var.ori, rho, rho.which,gene
      )%>%#,-rho.bw.tp,-p.bwtp) %>% 
    spread(key=dtname,value=val) %>%
    rename(var=var.ori,x=trc,y=prt),
  .(id),
  scat.plot_permu.simu,'mRNA','Protein') # ALL GENES
dev.off()

pdf(
  sprintf(
    '%s%s_%s_%s_bwtp.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    scat_plot.pdf.prefix
  )
  ,width=16,
  height=10
)
scat.plot.list <- dlply(
  sim.dat.sample_rho.est %>% 
    dplyr::select(
      val, id, dtname, var.ori, rho, rho.which,gene
      ) %>%
    inner_join(
      sim.dat.sample_rho.est %>% 
        dplyr::select(
          val, id, dtname, var.new, rho, rho.which,gene
          ) ,
      by=c(setNames('var.new','var.ori'),'id','dtname','gene','rho.which','rho')
    ) %>%
    dplyr::rename( x=val.x, y=val.y, var=var.ori) %>%
    mutate(gene=paste(gene,dtname,sep='_')),
  .(rho.which,id), 
  scat.plot_permu.simu,
  x.lab='TimePoint = t-1',y.lab='TimePoint = t') # ALL GENES
dev.off()

#```
#Empirical test of Rho /n

#NULL distribution of each cluster is created /n
#by the frequency of Rho from `r itt_of_permute` permutations)

#```{r empirical test for Rho using permutation test  }


#//////////////////
# Permutation Test
#//////////////////

empi.rho.2 <- function(data,startCol){
  # Input data is required following variables ; 
  #  id, dtname, var, val, 
  #
  for(i in 1:itt_of_permute){
    if(i!=1){
      mstSfl <- ddply(    
        ## apply 'varSfl'(Shufling variables) function 
        ## to data by 'id'*'dtname'
        data,.(id,dtname),
        varSfl,sfl=TRUE,startCol=startCol,endCol=length(data))#,.progress='text')
        ## data must be numeric in cols from 'startCol' to 'endCol' 
    } else {
      mstSfl <- data
      ## if first itteration by 'id' or 'dtname' then don't varSfl
    }
    if (rho.interTemporal==TRUE){
      mstSflLong <- mstSfl %>%
        gather(var,val,starts_with('d_')) %>%
        spread(dtname,val) %>%
        mutate(
          var.prt = match( var, 
                           names(
                             data %>%
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
                             data %>%
                               dplyr::select(
                                 one_of(cols)
                               )
                           )
          ),
          var.trc = var.prt # + delay ## Comment-out modification for Simu.PG
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


sim.dat <- rbind(
  rbind(
    sim_norm_10000_trc,
    sim_norm_10000_prt
  )
  ) %>%
  mutate(id = paste('SimNo.',count2,sep='')) %>%
  dplyr::select(-count2) 


output.empi.rho <- ddply(sim.dat %>%
                           dplyr::select(
                             -var.new,-var.ori#,
                             #-rho.bw.tp,-p.bwtp
                             ) %>% 
                           spread(key=var,value=val) %>%
                           rename(
                             edges=id,
                             id=gene
                             ) %>%
                           dplyr::select(
                           #rho.which,
                           edges,id,dtname,one_of(cols)
                           ),
                         .(
                         #rho.which,
                         edges),
                         empi.rho.2, startCol=5,
                         .progress = 'text'
                         ) %>%
  ddply(
          .(
           #rho.which ,
           edges
           ) ,
          my.rank ,
          'rho.est'
          ) %>%
  mutate(
    p_value = 1-(rank / itt_of_permute)
  ) 

Q.value_for_Rho <- function(output.empi.rho){
  output.empi.rho.which.q <- output.empi.rho %>%
    filter(itt==1 ) %>%
     # group_by(
     # rho.which
     #)
     # %>%
    mutate(
      q_value          = p.adjust(p_value,method = 'BH'),
      q_value.spearman = p.adjust(pvalue ,method = 'BH')
    ) %>%
    mutate(
      Flg_FDR.pt = ifelse(
        q_value < 0.05,1,0
      ),
      Flg_FDR.sp = ifelse(
        q_value.spearman < 0.05,1,0
      ),
      Flg_pval.pt = ifelse(
        p_value < 0.05,1,0
      ),
      Flg_pval.sp = ifelse(
        pvalue < 0.05,1,0
      )
    )  %>%
    ungroup() %>%
    data.frame()
  return(output.empi.rho.which.q)
}
output.empi.rho.which.q <- ddply(
  output.empi.rho,
  .(
   #rho.which
   ),
  Q.value_for_Rho
  )

output.empi.rho.which.q <- inner_join(
    output.empi.rho,
    output.empi.rho.which.q %>%
      dplyr::select(
        edges,
        #rho.which,
        q_value,
        q_value.spearman,
        starts_with('Flg')
        ),
    by=c(
     'edges'#,'rho.which'
     )
   )

#---

write.csv2(
  output.empi.rho.which.q ,
          file=sprintf(
            '%s%s_%s_%s_Q.val.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            permuted_Rho.csv.prefix
          ) 
)
write.csv2(
  output.empi.rho.which.q %>%
    filter(itt==1),
  file=sprintf(
    '%s%s_%s_%s_Q.val_observed.csv',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    permuted_Rho.csv.prefix
  ) 
)

#```

#```{r Histogram}
## Histogram

qtl.empi.rho <- output.empi.rho.which.q %>%
  ddply(
   .(
    edges #,
    #rho.which
    ),
    my.qtl,
    'rho.est'
    ) %>%
  rename(
    Q.025=V1,
    Q.050=V2,
    Q.950=V3,
    Q.975=V4
  ) %>%
  gather(var,val,Q.950)#starts_with('Q'))

df.perm_null_rho_histo <- inner_join(
  qtl.empi.rho,
  output.empi.rho.which.q %>%
    distinct(
     #rho.which,
     edges,
     itt
     ) %>%
    dplyr::select(
      starts_with('rho'), 
      edges,
      itt,
      starts_with('p'),
      starts_with('q'),
      starts_with('Flg')
      ),
  by=c(
   #'rho.which',
   'edges')
  ) #%>%
#  filter(
#    Flg_pval.pt==1 | Flg_pval.sp==1
#  )

perm_null_rho_histo <- function( data){
  
  edges <- unique(data$edges)
  
  p     <- ggplot(data ,aes(x=rho.est))
  p     <- p + geom_histogram(position='identity',binwidth = 0.01)
  
  gVline   <- geom_vline(
    aes(
      xintercept = rho.est,
      colour     = 'red',
      linetype='Observed'),
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
  
  
  legend.title <- labs(
    title=paste('null dist. of Rhos from ',
                itt_of_permute,
                ' permutations;',edges,sep=''),
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
    '%s%s_%s_%s.pdf',
    outputDirectry.prefix,
    outprefix,
    bioproc,
    perm_null_rho_histo.pdf.prefix
  )
)
perm.rho.histo <- dlply(
  df.perm_null_rho_histo,
  .(
   #rho.which,
   edges
   ),
  perm_null_rho_histo
)
dev.off() ## pdf(perm_null_rho_histo)

#```
