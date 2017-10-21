
# Dependency: func_for_calcEscore_20170225 or later version

dlply(dataGQ,.(id),gq_method2)

## the pentose phosphate pathway
##
# 5-pospho-ribosyl-1(alpha)-pyrophosphate synthetase
gq_method(dataGQ %>%
            filter(id=='YBL068W', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YBL068W', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)

# Ribulose-phosphate 3-epimerase
gq_method(dataGQ %>%
            filter(id=='YJL121C', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YJL121C', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)

# Protein with a possible role in t-RNA export
gq_method(dataGQ %>%
            filter(id=='YNR034W', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YNR034W', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)


## the glycine serine threonine metabolism pathway
##

# T subunit of the mitochondrial glycine decarboxylase complex
gq_method(dataGQ %>%
            filter(id=='YDR019C', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YDR019C', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)

# ATP-binding protein of unknown function
gq_method(dataGQ %>%
            filter(id=='YGR205W', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YGR205W', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)

# 3-phosphoglycerate dehydrogenase and alpha-ketoglutarate reductase
gq_method(dataGQ %>%
            filter(id=='YIL074C', dtname=='trc') %>%
            dplyr::select(starts_with('d_')),
          dataGQ %>%
            filter(id=='YIL074C', dtname=='prt') %>%
            dplyr::select(starts_with('d_')),
          timepoint = 4)
