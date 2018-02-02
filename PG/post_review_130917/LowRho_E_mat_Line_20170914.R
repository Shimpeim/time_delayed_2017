source('/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Functions/func_data_mnp_20151125.R')

data.path <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/Data/Selevsek_2015_all_Selevsek_20160706.RData'
data <- load(data.path)

source('/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Functions/require_packages_29Aug17.R')
Bibtex.out(Bibtex=TRUE)

source('/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/Functions/func_for_calcEscore_20170225.R')
source('./PG/Functions/func_for_permutest_of_Rho_20160819.R')

output.path <- '/Users/mos/Dropbox/Draft_201603Morimoto/Analysis/PG/post_review_130917'

z_norm        = TRUE
interTemporal = TRUE

if(z_norm == TRUE){
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

#==== E-matrix ====#

# Cluster ID#16
a_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YDR019C')) # =
_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YGR205W')) # =
_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YIL074C')) # =

# Cluster ID#16
VPS54_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YDR027C')) # = VPS54
RPL26A_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YLR344W')) # = RPL26A


# Cluster ID#16
VPS54_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YDR027C')) # = VPS54
RPL26A_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YLR344W')) # = RPL26A

# Cluster ID#15
SRM1_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YGL097W')) # = SRM1
UTP14_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YML093W')) # = UTP14

# Cluster ID#14
RLP24_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YLR009W')) # = RLP24
DCS1_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YLR270W')) # = DCS1

# Cluster ID#9
GND2_e  <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YGR256W')) # = GND2
RPL6A_e <- gq_scoreMatOut_withName(dataGQ %>% filter(id%in%'YML073C')) # = RPL6A

write.csv(VPS54_e[[1]],file=sprintf('%s/%s.csv',output.path,'VPS54_escore'))
write.csv(RPL26A_e[[1]],file=sprintf('%s/%s.csv',output.path,'RPL26A_escore'))
write.csv(SRM1_e[[1]],file=sprintf('%s/%s.csv',output.path,'SRM1_escore'))
write.csv(UTP14_e[[1]],file=sprintf('%s/%s.csv',output.path,'UTP14_escore'))
write.csv(RLP24_e[[1]],file=sprintf('%s/%s.csv',output.path,'RLP24_escore'))
write.csv(DCS1_e[[1]],file=sprintf('%s/%s.csv',output.path,'DCS1_escore'))
write.csv(GND2_e[[1]],file=sprintf('%s/%s.csv',output.path,'GND2_escore'))
write.csv(RPL6A_e[[1]],file=sprintf('%s/%s.csv',output.path,'RPL6A_escore'))

#==== Line-plot ====#

# levels for x-axis #

cols <- c('0min','30min','60min','90min','120min')
data_ana_norm_for_line <- data_ana_norm %>%
  mutate(
    var=factor(var,cols)
  )

# function #

line.plot_lite <- function(data){
  id <- data$id[1]
  p <- ggplot(data,
              aes(x=var,y=val,group=dtname)
  )
  plot.type <- geom_line(aes(colour=dtname),size=1.5)
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
        'Gene name_',id,sep=''
        )
      )
    )
  plot.ylab <- ylab('Normalised')
  plot( p+ plot.type+ plot.theme + plot.title + plot.ylim + plot.ylab)
}

pdf(
  sprintf(
    '%s/%s.pdf',
    output.path,
    'LowRho_line_plot'
  ),
  width=16,
  height=10
)

line.plot.list <- 
  dlply(
    data_ana_norm_for_line %>% 
      filter(
        id %in% c(
          # Cluster 16
          'YDR027C', #
          'YLR344W', #
          # Cluster 15
          'YGL097W', #
          'YML093W', #
          # Cluster 14
          'YLR009W', #
          'YLR270W', #
          # Cluster 9
          'YGR256W', #GND2
          'YML073C'  #RPL6A
          ),
        var %in% cols
        ),
    .(id),
    line.plot_lite
  ) 
dev.off()
