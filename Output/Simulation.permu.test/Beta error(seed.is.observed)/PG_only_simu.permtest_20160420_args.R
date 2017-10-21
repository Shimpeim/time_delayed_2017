#R --vanilla --quiet < PG_only_simu.permtest_20160411_args.R --args Selevsek_2015_all_Selevsek_20160123.RData  YDL110C_YDR122W 100 0.3 0.05  > Simulation.permtest_YDL110C_YDR122W.log 2>&1 &
#R --vanilla --quiet < PG_only_simu.permtest_20160411_args.R --args Selevsek_2015_all_Selevsek_20160123.RData  YBR230C_YHL021C 100 0.3 0.05 > Simulation.permtest_YBR230C_YHL021C.log 2>&1 &
#R --vanilla --quiet < PG_only_simu.permtest_20160411_args.R --args Selevsek_2015_all_Selevsek_20160123.RData  YBR035C_YHR179W_YJL217W 100 0.3 0.05 > Simulation.permtest_YBR035C_YHR179W_YJL217W.log 2>&1 &
#R --vanilla --quiet < PG_only_simu.permtest_20160411_args.R --args Selevsek_2015_all_Selevsek_20160123.RData  YBR016W_YMR175W 100 0.3 0.05 > Simulation.permtest_YBR016W_YMR175W.log 2>&1 &
#R --vanilla --quiet < PG_only_simu.permtest_20160411_args.R --args Selevsek_2015_all_Selevsek_20160123.RData  YBR016W_YMR175W 100 0.3 0.05 > Simulation.permtest_YBR016W_YMR175W.log 2>&1 &

# Cluster-wise simulation for the testing type-I & II error of the permutation test
# by adding random noise to the time-series data of genes' expression quantity log2FC.   
# Morimoto,S. 2016

#This source code is built under R ver.3.1.2

c_args <- commandArgs(trailingOnly=T)
makedataData  <- c_args[1]
c_id.select   <- c_args[2]
c_sampleSize  <- c_args[3]
c_rdm.err     <- c_args[4]
c_rho.intval  <- c_args[5]

sampleSize <- as.numeric(c_sampleSize) # N of random sampled data for each rho 
rdm.err    <- as.numeric(c_rdm.err)    # Define the SD of the random error 
rho.intval <- as.numeric(c_rho.intval) # Define the bin around each rho 

rho.seq        <- c(seq(0.0,0.3,by=0.1))         #
rho.seq.bw.tp  <- c(seq(0.0,0.3,by=0.1))         #

outprefix      <- sprintf('%s_sd_%s','simulation_permtest',c_rdm.err)
bioproc        <- c_id.select
outputDirectry.prefix          <- './'
#outputDirectry.prefix         <- '/Users/mos/Dropbox/D/Thesis001/PG/Simulation/Permu.test/'

dataDirectry <- './'
funcDirectry <- './'
pkgsDirectry  <- './'


#makedataData   <- 'Selevsek_2015_all_Selevsek_20160123.RData'

#dataDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test/'
#funcDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test'
#pkgsDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Simulation/Permu.test'


rho.interTemporal <- FALSE # convert data for corr.Rho calc. to inter-temporal ?

#```{r prefix for output files}



histogram.pdf.prefix        <- 'histogram'

scat_plot.pdf.prefix        <- 'scatt_'
line_plot.pdf.prefix        <- 'line_'

perm_null_rho_histo.pdf.prefix <- 'permNullDistHisto'
permuted_Rho.csv.prefix        <- 'permuted_Rho'

RData_PermRhoAnalysis_save.image.prefix <- 'permuted_Rho'


# permuted Rho

itt_of_permute <- 10000
nboot <- 100

#```

#```{r functions sorce code file}

permtRho   <- 'func_for_permutest_of_Rho_20151124.R' 
dataManu   <- 'OldFunc_20151125.R'
EscoreCalc <- 'func_for_calcEscore_20160123.R' # modified: 2016/01/23
#```


#```{r Load libraries}

## LIBRARIES

packages <- c(
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


#```{r load R_data}

load(
  file = sprintf(fmt = '%s%s',dataDirectry,makedataData)
)
sim.dat.seed <- data_ana %>%
  filter(id %in% id.select)

w.timePoint  <- timePoint

cols <- unique(data_ana$var)
bioproc        <- c_id.select
#```

#```{r generating sim.dat}

##========================##
##  Rho of observed data  ##
##========================##

sim.dat.seed_1 <- sim.dat.seed %>%
  group_by(id,dtname) %>%
  mutate(
    norm = scale(val, center = TRUE,scale = TRUE),
    var.ori=match(var,cols),
    var.new=var.ori+1
  ) %>%
  ungroup()

## Rho between y_mRNA and y_protein ( of Observed Data ) ##

rho.ori <- sim.dat.seed_1 %>%
  mutate(keykey = paste(id,var,sep='_')) %>%
  dplyr::select(-id,-var,-val,-var.ori,-var.new) %>% 
  spread(key = keykey,value=norm) %>%
  dplyr::select(-dtname) %>%
  t() %>% 
  as.matrix() %>%
  cor(method='s')

## Rho between y_(t) and y_(t-1) ( of Observed Data ) ##


rho.ori_inter.tmprl <- sim.dat.seed_1 %>%
  select(id,dtname,norm,var.ori)%>%
  inner_join(
    sim.dat.seed_1 %>%
      select(id,dtname,norm,var.new) %>%
      rename(var.ori=var.new),
    by=c('id','dtname','var.ori')
    ) %>%
  dplyr::select(norm.x,norm.y) %>%
  as.matrix() %>%
  cor(method='s')
##====================================================##

##====================##
## Simulation dataset ##
##====================##

count     <- 0 # Initiate the counter in 'repeat{do if(count==sampleSize){break}}'
repeat{  
  
  # Do until the condition 'count == sampleSize '
  # (as written in if statement) is met.
  
  sim.dat <- sim.dat.seed %>%
    
    rename(w.val=val) %>%
    mutate(
      dev = rnorm(length(id),0,rdm.err),  # New variable 'dev' 
      # is created from rnorm 
      val = ifelse(var=='T1-T1',w.val,w.val+dev),
      # This conditioning is for data 
      # such that the timepoint:0 (Control) is always 0
      var.ori=match(var,cols),
      var.new=var.ori+1
    ) %>%
    dplyr::select(id,dtname,var,val,var.ori,var.new) %>%
    group_by(id,dtname) %>%               # statement the group for Z-standardization 
    mutate(norm = scale(val,center = TRUE,scale=TRUE)) %>% # Z-standardization
    ungroup()                             # Ungroup 'group_by(id,dtname)'
  
  ## Rho between y_mRNA and y_protein ( of Simulation Data ) ##
  
  sim.dat.rho <- sim.dat %>%
    # of the normalized data
    mutate(keykey = paste(id,var,sep='_')) %>%
    dplyr::select(-id,-var,-val,-var.ori,-var.new) %>% 
    spread(key = keykey,value=norm) %>%   # variable 'norm' is normalised 'val'
    dplyr::select(-dtname) %>%
    t() %>% 
    as.matrix() %>%
    cor(method='s')
  
  ## Rho between y(t) and y(t-1) ( of Simulation Data ) ##
  
  sim.dat.rho_inter.tmprl <- sim.dat %>%
    select(id,dtname,norm,var.ori)%>%
    inner_join(
      sim.dat %>%
        select(id,dtname,norm,var.new) %>%
        rename(var.ori=var.new),
      by=c('id','dtname','var.ori')
    ) %>%
    dplyr::select(norm.x,norm.y) %>%
    as.matrix() %>%
    cor(method='s')
  ## 
  
  rho.which <- which(
    abs(sim.dat.rho[1,2] + rho.seq - rho.ori[1,2] ) <  rho.intval
  )
  rho.which.bw.tp <- which(
    abs(sim.dat.rho_inter.tmprl[1,2] + rho.seq.bw.tp - rho.ori_inter.tmprl[1,2] 
        ) <  rho.intval
  )
  #print(abs(sim.dat.rho[1,2] + rho.seq - rho.ori[1,2] ))
  
  if(length(rho.which) > 0 & length(rho.which.bw.tp) > 0 ){
    
    sim.dat$rho.which       <- rho.which
    sim.dat$rho.which.bw.tp <- rho.which.bw.tp
    sim.dat$rho_combi <- paste(rho.which,rho.which.bw.tp,sep='_')  
    print(sim.dat$rho_combi)  
    sim.dat$rho       <- sim.dat.rho[1,2]
    count <- count + 1  # counter of
    if(count==1){
      sim.dat$count   <- count
      sim.dat.set     <- sim.dat
    }else{
      sim.dat$count <- count
      sim.dat.set   <- rbind(sim.dat.set, sim.dat)
    }


    if(min(plyr::count(sim.dat.set,'rho_combi')[,2]
           ) == sampleSize*2*timePoint*length(id.select)){
      break
    }
  }
}

sim.dat.set.smpler <- sim.dat.set %>%
  dplyr::select(rho.which,count) %>%
  distinct(count)%>%
  ddply(
    .(rho.which),
    sample_n,
    sampleSize,
    replace = FALSE
    ) 
write.csv(sim.dat.set.smpler,
          file = sprintf(
            '%s%s_%s_%s.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            '_sampler')
)
sim.dat.set.smpl <- sim.dat.set %>%
  inner_join(
    sim.dat.set.smpler,
    by=c('count','rho.which')
    )

write.csv(sim.dat.set.smpl,
      file = sprintf(
        '%s%s_%s_%s.csv',
        outputDirectry.prefix,
        outprefix,
        bioproc,
        'simu_dat')
)


#```{r make pdf file of time series line plots}
line.plot <- function(data){
  id <- data$id[1]
  edge  <- data$count[1]
  rho.which <- data$rho.which
  rho       <- round(data$rho,2)
  p <- ggplot(data %>%
                mutate(timepoint=paste('timepoint ',var)),
              aes(x=var,y=norm,group=dtname)
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
        'Gene name_',id,'[rho.which:',rho.which,'(Rho=',rho,'),count:',edge,']',
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
    sim.dat.set.smpl ,
    .(rho_combi,count,id),
    line.plot
  ) 
dev.off()


#```{r plotting scatter plot}

write.csv(
  sim.dat.set.smpl %>%
    dplyr::select(-val) %>%
    spread(key=var,value=norm),
          file = sprintf(
            '%s%s_%s_%s.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            'spread')
  )


scat.plot_permu.simu <- function(data){
  shape_basket <- c(8,15,16,17,18,19) # plot shape selecter. (ref. help(shape) {ggplot2})
  edge  <- data$count[1]
  RHO   <- cor.test_by.id(data,'prt','trc')$rho.est
  rho.which <- data$rho.which[1]
  p <- ggplot(data %>%
                mutate(timepoint=paste('timepoint #',var),
                       cluster=paste('Group;',count)
                ),
              aes(x=trc,y=prt,group=timepoint))
  
  if(length(unique(data$id))<=6){
    plot.type <- geom_point(aes(colour=timepoint,shape=id),size=10,alpha=0.6)
  }else{
    plot.type <- geom_point(aes(colour=timepoint,shape=id),size=4,alpha=0.6)
  }
  
  if(length(unique(data$id))<=6){
    plot.shape <- scale_shape_manual(values=shape_basket[c(1:length(unique(data$id)))])
  }else{
    plot.shape <- scale_shape_manual(values=rep(18,length(unique(data$id))))
  }
  
  legends.guide <-guides(
    guide_legend(keywidth = 0.3, keyheight = 0.3,nrow = 2, byrow = TRUE)
  )
  plot.title <- labs(
    list(
      title = paste(
        '[rho.which:',rho.which,'(Rho=',RHO,'),count:',edge,']',
        sep=''
      )
    )
  )
  plot.bg <- theme_bw()
  plot( p+ plot.type+ plot.shape + plot.title + plot.bg +legends.guide)
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
  sim.dat.set.smpl %>% 
    dplyr::select(-val) %>% 
    spread(key=dtname,value=norm),
  .(rho.which,count), 
  scat.plot_permu.simu) # ALL GENES

dev.off()

#```
#Empirical test of Rho /n

#NULL distribution of each cluster is created /n
#by the frequency of Rho from `r itt_of_permute` permutations)

#```{r empirical test for Rho using permutation test  }

# make data for permutest

empi.rho <- function(data){
  for(i in 1:itt_of_permute){
    if(i!=1){
      mstSfl <- ddply(
        data,.(id,dtname),
        varSfl,sfl=TRUE,startCol=3,endCol=length(data))#,.progress='text')
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


output.empi.rho <- ddply(sim.dat.set.smpl %>%
                           dplyr::select(-val) %>% 
                           spread(key=var,value=norm),
                         .(count),empi.rho,.progress = 'text') %>%
  ddply(.(count),my.rank,'rho.est') %>%
  mutate(
    p_value = 1-(rank / itt_of_permute)
  ) 

panderOptions('table.split.table', 1000)
panderOptions('table.continues', '')
pander(
  output.empi.rho %>%
    filter(itt==1 ) %>%
    dplyr::select(rho.est,p_value,count)#,V1,V4)
)

output.empi.rho.q <- output.empi.rho %>%
  filter(itt==1 ) %>%
  inner_join(sim.dat.set.smpler,by=c('count')
             ) #%>%
#  mutate(q_value=p.adjust(p_value,method = 'BH'))

write.csv(output.empi.rho.q,
          file=sprintf(
            '%s%s_%s_%s.csv',
            outputDirectry.prefix,
            outprefix,
            bioproc,
            permuted_Rho.csv.prefix
          ) 
)

#```

#```{r Histogram}
## Histogram

head(output.empi.rho)
qtl.empi.rho <- output.empi.rho %>%
  ddply(.(count),my.qtl,'rho.est') %>%
  rename(
    Q.025=V1,
    Q.050=V2,
    Q.950=V3,
    Q.975=V4
  ) %>%
  gather(var,val,Q.950)#starts_with('Q'))

df.perm_null_rho_histo <- inner_join(qtl.empi.rho,output.empi.rho)

perm_null_rho_histo <- function( data){
  
  count <- unique(data$count)
  
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
  
  
  legend.title <- labs(
    title=paste('null dist. of Rhos from ',
                itt_of_permute,
                ' permutations; count=',count,sep=''),
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
  .(count),
  perm_null_rho_histo
)
dev.off() ## pdf(perm_null_rho_histo)

#```
