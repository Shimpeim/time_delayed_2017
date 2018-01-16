#```{r BioMart setting}

# Gene id convert (via bioMart)
# 
#bioMartDB      <- "fungi_mart_29" # 2015.10.28 DB update                                     # biomaRt::listMarts()
#bioMartDataSet <- 'scerevisiae_eg_gene'                     
#                                  # biomaRt::listDatasets(db)
#inputName      <- 'wikigene_name'   # 'wikigene_name' or 'ensembl_gene_id' 
#                                  # biomaRt::listFilters(sceg)
#outputName     <- 'uniprot_swissprot_accession'
#                                  # biomaRt::listFilters(sceg)

bioMartHost    <- 'fungi.ensembl.org' # 2015.11.10 updated (HOST='biomart.org' has stopped)
bioMartDB      <- "fungal_mart" # 2015.11.10 DB update    
# biomaRt::listMarts(host=bioMartHost)
bioMartDataSet <- 'scerevisiae_eg_gene'                     
# biomaRt::listDatasets(db)
inputName_1      <-  'ensembl_gene_id' # 'wikigene_name' or 'ensembl_gene_id' 
# biomaRt::listFilters(sceg)
outputName     <- 'uniprot_swissprot_accession'
# biomaRt::listFilters(sceg)s

#```