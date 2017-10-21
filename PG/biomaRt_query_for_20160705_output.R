packages <- c(
  'dplyr', # progress bar
  'plyr',  # progress bar
  'dplyr', # progress bar
  'tidyr',
  'xlsx',
  'readxl',
  'ggplot2',
  'gplots',
  'GMD',
  'pvclust',
  'reshape2',
  'pander',
  'stringr'
#  'biomaRt'
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
require('xlsx')
require('readxl')
require('ggplot2')
require('gplots')
require('GMD')
require('pvclust')
require('reshape2')
require('pander')
require('stringr')if(  !('biomaRt' %in% installed.packages())){  source('http://bioconductor.org/biocLite.R')  biocLite('biomaRt')   require('biomaRt')}else   require('biomaRt')
  dataDirectry <- '/Users/mos/Dropbox/Draft_201603Morimoto/Work_Fig'dataPath_1     <- 'HighAUcluster_20160705.csv'bioMartHost    <- 'fungi.ensembl.org'# 2015.11.10 updated (HOST='biomart.org' has stopped)bioMartDB      <- "fungal_mart" # 2015.11.10 DB update    # biomaRt::listMarts(host=bioMartHost)bioMartDataSet <- 'scerevisiae_eg_gene'                     # biomaRt::listDatasets(db)inputName_1    <-  'ensembl_gene_id' # 'wikigene_name' or 'ensembl_gene_id' # biomaRt::listFilters(sceg)outputName     <- 'description'# biomaRt::listAttributes(sceg)db.DL <- useMart(bioMartDB,host=bioMartHost) # listMarts(host=bioMartHost)sceg <- useDataset(bioMartDataSet, mart = db.DL) # listDatasets(db)

HighAUcluster <- read.csv(sprintf('%s/%s', dataDirectry,dataPath_1),header = TRUE)

merged_1 <- HighAUcluster %>%
    inner_join(
        getBM(
            attributes = c(inputName_1,outputName,'external_gene_id'),
            filters =  c(inputName_1), # listFilters(sceg)
            values = HighAUcluster$V1,#HighAUcluster$V2),#HighAUcluster$id,
            mart = sceg
            ),
        by=c('V1'='ensembl_gene_id')
    ) %>%
    inner_join(
        getBM(
            attributes = c(inputName_1,outputName,'external_gene_id'),
            filters =  c(inputName_1), # listFilters(sceg)
            values = HighAUcluster$V2,#HighAUcluster$V2),#HighAUcluster$id,
            mart = sceg
            ),
        by=c('V2'='ensembl_gene_id')
    )%>%
          mutate(
              description.xy = paste(
                paste( external_gene_id.x,description.x, sep=' : \r'),
                paste( external_gene_id.y,description.y, sep=' : \r'),
                '\r',
                sep=c('\r\r')
                  )
             ) 


write.csv(merged_1,
          file=sprintf(fmt = '%s/%s.csv',
                  dataDirectry,
                  'genename_protein_description')
          ) 



id_convert[,'id']    <- id_convert[,inputName_1]
id_convert[,'description'] <- id_convert[,outputName]

write.csv(id_convert %>%
            dplyr::select(id,uniid),
          file=sprintf(fmt = '%s/%s.csv',
                  dataDirectry,
                  'protein_description')
          )