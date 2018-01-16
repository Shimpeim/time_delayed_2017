#---
#  title: "Application Qian 2001 method on Cheng-Vogel 2016 data (high corr b/w rep 1 and 2)"
#author: "s_morimoto"
#date: "2017/12/27"
#output: 
#  pdf_document: 
#  latex_engine: xelatex
#---
#  This source code is built under R ver.3.1.2

#You are using  `r R.version.string` 



#SETUP parameters
#```{r setting up parameters}
rm()

##--- setting : data  ##

gene_filter <- "corr_gt_0.7_Pier ==1"

dataDirectry   <- "./Cheng_Vogel/Data"

dataPath_trc     <- "rep1_rep2_corr.xlsx"
#'Fournier_20151018_trc.xlsx'
#"Gene2011_testResult_ana.xlsx" # 
dataPath_prt     <- "rep1_rep2_corr.xlsx"
#'Fournier_20151018_prt.xlsx'
#"ProtSWATH_testResult_ana.xlsx" # 
#dataPath_pathway <- "2015_Selevsek_proteins_KEGG.xlsx" #  
bioproc <- NULL #'pentose' 'GlySerThrMetabo' #'pentose'# #NULL#

typeBoolean  <- FALSE
sheet    = 2
#sheetPathway = 1

header   = TRUE

startRow_trc = 1
startRow_prt = 1
#startRow_pw  = 1

endRow_trc = 1238 
endRow_prt = 1238 
#endRow_pw  = 54#54

timePoint <- 8
repet     <- 1

discrete_x_label <- c('normalized relative expression')#'inter-temporal')
continuous_y_label <- c('0h is 1')


cols <- c('t_0h','t_0.5h','t_1h','t_2h','t_8h','t_16h','t_24h','t_30h')

## end of setting : data  ---##

columnName <- c(rep(cols,repet))

## --- setting : data processing before analysis #

z_norm <- TRUE  ## Z-Normalize by genes and datanames ( TRUE/ FALSE)
interTemporal <- FALSE

## end of setting : data processing before analysis ---##

## --- setting : filenames of the OUTPUT FILEs ##

allOutput.prefix <- 'Cheng_Vogel_2016_rep_mean'

outputDirectry <- 
  "./Cheng_Vogel/Output"

file_dataGQ.prefix         <- 'dataGQ_Cheng_Vogel'
file_data_idConvert.prefix <- 'id.convert_Cheng_Vogel'
file_save.image     <- 'Cheng_Vogel_rep_PC07_.RData'

## end of setting : filenames of the OUTPUT FILEs ---##

#```

#LOAD PACKAGES
#```{r install packages,echo=FALSE}

funcDirectry <- './PG/Functions'
source(file=sprintf('%s/%s',funcDirectry,'require_packages_25NOV17.R'))

## avoid the error (source: https://code.google.com/p/rexcel/issues/detail?id=33)
#options(java.parameters = "-Xmx1000m") 

panderOptions('table.split.table', 1000)
panderOptions('table.continues', '')

#```

#LOAD DATA
#```{r loading data }

raw_data.trc <- read_excel(
  sprintf(fmt = '%s/%s',
          dataDirectry,
          dataPath_trc
  ),
  sheet,
  col_names = header,
  col_types = NULL, 
  na = "",
  skip = startRow_trc-1
  ) %>%
  filter(
    eval(
      parse(
        text=gene_filter
        )
      )
    ) %>%
  dplyr::select(
    'Gene', 'GeneSymbol', 'ProteinNames', 'Cluster', 
    "R1t0.5_1",
    "R1t1_1",
    "R1t2_1",
    "R1t8_1",
    "R1t16_1",
    "R1t24_1",
    "R1t30_1"
    ) %>%
  mutate(R1t0_1=0) %>%
  rename(
    t_0h   = R1t0_1,
    t_0.5h = R1t0.5_1,
    t_1h   = R1t1_1,
    t_2h   = R1t2_1,
    t_8h   = R1t8_1,
    t_16h  = R1t16_1,
    t_24h  = R1t24_1,
    t_30h  = R1t30_1
  ) %>%
  gather(var,log_val,-Gene, -GeneSymbol, -ProteinNames, -Cluster) %>%
  mutate(
    val_rep1 = 2**log_val, 
    var=factor(var,levels=cols)
    ) %>%
  rename(
    id=Gene) %>%
  inner_join(
    read_excel(
      sprintf(fmt = '%s/%s',
              dataDirectry,
              dataPath_trc
      ),
      sheet,
      col_names = header,
      col_types = NULL, 
      na = "",
      skip = startRow_trc-1
    ) %>%
      filter(
        eval(
          parse(
            text=gene_filter
          )
        )
      ) %>%
      dplyr::select(
        'Gene', 'GeneSymbol', 'ProteinNames', 'Cluster', 
        "R2t0.5_1",
        "R2t1_1",
        "R2t2_1",
        "R2t8_1",
        "R2t16_1",
        "R2t24_1",
        "R2t30_1"
      ) %>%
      mutate(R2t0_1=0) %>%
      rename(
        t_0h   = R2t0_1,
        t_0.5h = R2t0.5_1,
        t_1h   = R2t1_1,
        t_2h   = R2t2_1,
        t_8h   = R2t8_1,
        t_16h  = R2t16_1,
        t_24h  = R2t24_1,
        t_30h  = R2t30_1
      ) %>%
      gather(var,log_val,-Gene, -GeneSymbol, -ProteinNames, -Cluster) %>%
      mutate(
        val_rep2 = 2**log_val,
        var=factor(var,levels=cols)
        ) %>%
      rename(
        id=Gene),
    by = c('id', 'GeneSymbol', 'ProteinNames', 'Cluster', 'var')
  ) %>%
  mutate(
    val = 1/2 *(val_rep1+val_rep2)
  )


raw_data.prt <- read_excel(
  sprintf(fmt = '%s/%s',
          dataDirectry,
          dataPath_trc
  ),
  sheet,
  col_names = header,
  col_types = NULL, 
  na = "",
  skip = startRow_trc-1
  ) %>%
  filter(
    eval(
      parse(
        text=gene_filter
      )
    )
  ) %>%
  dplyr::select(
    'Gene', 'GeneSymbol', 'ProteinNames', 'Cluster',
    "LFQ.intensity.2_05h_RS1",
    "LFQ.intensity.3_1h_RS1",
    "LFQ.intensity.4_2h_RS1",
    "LFQ.intensity.5_8h_RS1",
    "LFQ.intensity.6_16h_RS1",
    "LFQ.intensity.7_24h_RS1",
    "LFQ.intensity.8_30h_RS1"
    ) %>%
  mutate(LFQ.intensity.1_0h_RS1=0) %>%
  rename(
    t_0h   = LFQ.intensity.1_0h_RS1,
    t_0.5h = LFQ.intensity.2_05h_RS1,
    t_1h   = LFQ.intensity.3_1h_RS1,
    t_2h   = LFQ.intensity.4_2h_RS1,
    t_8h   = LFQ.intensity.5_8h_RS1,
    t_16h  = LFQ.intensity.6_16h_RS1,
    t_24h  = LFQ.intensity.7_24h_RS1,
    t_30h  = LFQ.intensity.8_30h_RS1
  ) %>%
  gather(var,log_val,-Gene, -GeneSymbol, -ProteinNames, -Cluster) %>%
  mutate(
    val_rep1 = 2**log_val,  
    var=factor(var,levels=cols)
    ) %>%
  rename(id=Gene) %>%
  
  inner_join(
    read_excel(
      sprintf(fmt = '%s/%s',
              dataDirectry,
              dataPath_trc
      ),
      sheet,
      col_names = header,
      col_types = NULL, 
      na = "",
      skip = startRow_trc-1
    ) %>%
      filter(
        eval(
          parse(
            text=gene_filter
          )
        )
      ) %>%
      dplyr::select(
        'Gene', 'GeneSymbol', 'ProteinNames', 'Cluster',
        "LFQ.intensity.2_05h_RS2",
        "LFQ.intensity.3_1h_RS2",
        "LFQ.intensity.4_2h_RS2",
        "LFQ.intensity.5_8h_RS2",
        "LFQ.intensity.6_16h_RS2",
        "LFQ.intensity.7_24h_RS2",
        "LFQ.intensity.8_30h_RS2"
      ) %>%
      mutate(LFQ.intensity.1_0h_RS2=0) %>%
      rename(
        t_0h   = LFQ.intensity.1_0h_RS2,
        t_0.5h = LFQ.intensity.2_05h_RS2,
        t_1h   = LFQ.intensity.3_1h_RS2,
        t_2h   = LFQ.intensity.4_2h_RS2,
        t_8h   = LFQ.intensity.5_8h_RS2,
        t_16h  = LFQ.intensity.6_16h_RS2,
        t_24h  = LFQ.intensity.7_24h_RS2,
        t_30h  = LFQ.intensity.8_30h_RS2
      ) %>%
      gather(var,log_val,-Gene, -GeneSymbol, -ProteinNames, -Cluster) %>%
      mutate(
        val_rep2 = 2**log_val, 
        var=factor(var,levels=cols)
        ) %>%
      rename(id=Gene),
    by = c('id', 'GeneSymbol', 'ProteinNames', 'Cluster', 'var')
  ) %>%
  mutate(
    val = 1/2 *(val_rep1+val_rep2)
  )

#```

#```{r tranforming data}
data.trc <- mutate(raw_data.trc,
                   dtname='trc'
) %>%
  dplyr::select(id,var,val,dtname)
  

data.prt <- mutate(raw_data.prt,
                   dtname='prt'
) %>%
  dplyr::select(id,var,val,dtname)

if (z_norm==TRUE){
  dataLong <- bind_rows(data.trc,data.prt) %>%
    rename(val2=val) %>%
    group_by(id,dtname) %>%  # " dataLong <- as.data.frame(dataLong) "
    mutate(val=scale(val2, center = TRUE, scale = TRUE)) %>%
    dplyr::select(-val2)
}else 
  dataLong <- bind_rows(data.trc,data.prt) 

attributes(dataLong$val) <- NULL 
# attrs are created by scale function which causes errors when this data treated as data.frame
dataLong <- as.data.frame(dataLong) 
# ungroup the BY-groups created by  " %>% group_by(id,dtname)) "

#```

#Data summary
#```{r Data summary}

data_ana <- 
  dataLong

#data_ana$var <- factor(data_ana$var,levels=cols)

summ_data_ana <- data_ana %>%
  group_by(dtname,id) %>%
  summarise(
    n=n(),mean=mean(val),sd=sd(val),
    min=min(val),median=median(val),max=max(val))
pander(summ_data_ana)

#```

#```{r saving RData file}
if(is.null(bioproc)) bioproc <- 'all'

save(
  data_ana,
  #  dataGQ,
  timePoint,
  repet,
  bioproc,
  discrete_x_label,
  continuous_y_label,
  cols,
  #  gq_method,
  file = sprintf(fmt = '%s/%s_%s_%s',
                 outputDirectry,
                 allOutput.prefix,
                 bioproc,
                 file_save.image
  )
)

