
packages <- c(
  'dplyr', # progress bar
  'plyr',  # progress bar
  'dplyr', # progress bar
  'tidyr',
  'xlsx',
  'ggplot2',
  'gplots',
  'GMD',
  'pvclust',
  'reshape2',
  'pander',
  'matrixcalc'
  #  'scaleboot' 
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
require('ggplot2')
require('gplots')
require('GMD')
require('pvclust')
require('reshape2')
require('pander')
require('matrixcalc')

#if( !("biomaRt" %in% installed.packages()) ){
#  source('http://bioconductor.org/biocLite.R')
#  biocLite( 'biomaRt' )
#  require( 'biomaRt' )
#}else  require( 'biomaRt' )

if( !("scaleboot" %in% installed.packages()) ){
  install.packages(
    sprintf(fmt = '%s%s',
            pkgsDirectry,
            "scaleboot_0.3-3.tar.gz"  
            # ver. 16-May-2010 16:21
            # https://cran.r-project.org/src/contrib/Archive/scaleboot/
    ),
    repos = NULL, type = "source")
  require( 'scaleboot' )
}else require( 'scaleboot' )


# citation #

Bibtex.out <- function(Bibtex){
  if(Bibtex){
    write(toBibtex(citation()),file="CRAN")
    for(i in 1:length(packages)){
      write(toBibtex(citation(packages[i])),file=sprintf("Biblio/%s%s.bib",packages[i],"_CRAN"))
      }
    }
  }

#```


