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
require('stringr')

if(
  !('biomaRt' %in% installed.packages())){
  source('http://bioconductor.org/biocLite.R')
  biocLite('biomaRt') 
  require('biomaRt')}else 
  require('biomaRt')
  
dataDirectry   #Defined in main 
dataPath_1     #Defined in main 

source('./PG/Functions/_Cheng_Vogel_biomaRt_settingR.R')

## instruction for Ensembl gene mart: 
##   http://www.ensembl.org/info/data/biomart/biomart_r_package.html#biomartexamples


# db.DL <- useMart(bioMartDB,host=bioMartHost) # listMarts(host=bioMartHost)
db.DL <- useEnsembl(biomart="ensembl",GRCh=37)
sceg <- useDataset(bioMartDataSet, mart = db.DL) # listDatasets(db)

HighAUcluster <- read.csv(sprintf('%s/%s', dataDirectry, dataPath_1),header = TRUE)

merged_1 <- HighAUcluster %>%
    inner_join(
        getBM(
            attributes = c(inputName_1,outputName,'external_gene_name'),
            filters =  c(inputName_1), # listFilters(sceg)
            values = HighAUcluster$id,#HighAUcluster$V2),#HighAUcluster$id,
            mart = sceg
            ),
        by=c('id'='ensembl_gene_id')
    ) %>%
    inner_join(
        getBM(
            attributes = c(inputName_1,'description','external_gene_name'),
            filters =  c(inputName_1), # listFilters(sceg)
            values = HighAUcluster$id,#HighAUcluster$V2),#HighAUcluster$id,
            mart = sceg
            ),
        by=c('id'='ensembl_gene_id')
    )%>%
          mutate(
              description = paste(
                paste( external_gene_name.x,description, sep=': \r')
                  )
             ) 

detach(package:biomaRt)

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