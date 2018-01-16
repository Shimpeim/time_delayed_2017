rep_corr <- raw_data.trc %>%
  dplyr::select(
    corr_rank_RNA,
    corr_rank_prt
    )

hist(rep_corr$corr_rank_RNA)
hist(rep_corr$corr_rank_prt)


rep1_rep2_corr <- hist2d(
  x=rep_corr$corr_rank_RNA, 
  y=rep_corr$corr_rank_prt,
  nbins = 20,
  col=c("white", rainbow(50)
        )
  )

contour( 
  rep1_rep2_corr$x, 
  rep1_rep2_corr$y, 
  rep1_rep2_corr$counts, 
  nlevels=4 )
persp( 
  rep1_rep2_corr$x, 
  rep1_rep2_corr$y, 
  rep1_rep2_corr$counts,
  ticktype="detailed", 
  theta=300, 
  phi=30,
  expand=0.5, shade=0.5, col="cyan", ltheta=-10)

