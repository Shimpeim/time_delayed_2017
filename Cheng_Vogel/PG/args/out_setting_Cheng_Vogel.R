#R --vanilla --quiet < Morimoto_TDC.2016.main_20160819_args.R --args Selevsek_2015_all_Selevsek_20160706.RData             out.setting=F clustering=T ScaleBoot=F  lag.neg.ign=T  Selevsek_2015_id.convert_all.csv             > Selevsek_2015.all_Qian.log 2>&1 &


makedataData     <- 'Cheng_Vogel_2016_rep_mean_all_Cheng_Vogel_rep_PC07_.RData'
c_clustering     <- 'clustering=T'
c_ScaleBoot      <- 'ScaleBoot=F'
c_lag.neg.ign    <- 'lag.neg.ign=T'
id.convertData   <- 'genename_protein_description.csv'

outputDirectry.prefix <- './Cheng_Vogel/Output/'
dataDirectry <- './Cheng_Vogel/Data/'
funcDirectry <- './PG/Functions/'
pkgsDirectry  <- './PG/Packages/'