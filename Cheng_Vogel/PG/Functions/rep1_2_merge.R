require(tidyr)
require(tplyr)
require(dplyr)


names(output.empi.rho)

outputDirectry.prefix <- './Cheng_Vogel/Output/'
dataname_1 <- 'pier_corr_0.8/_Analyses_results_Cheng_Vogel_2016_rep1_all_Cheng_Vogel_rep_PC08__all_high_BP_uniProtAccess_Qian,2001_nboot10000_sb.off_.csv'
dataname_2 <- 'pier_corr_0.8/_Analyses_results_Cheng_Vogel_2016_rep2_all_Cheng_Vogel_rep_PC08__all_high_BP_uniProtAccess_Qian,2001_nboot10000_sb.off_.csv'

test <- read.csv(
  sprintf(
    '%s%s',
    outputDirectry.prefix,
    dataname_1
    )
  ) %>%
  full_join(
    read.csv(
      sprintf(
        '%s%s',
        outputDirectry.prefix,
        dataname_2
        )
      ),
    by = 'id'
  )

write.csv(
  test,
  sprintf(
    '%s%s',
    outputDirectry.prefix,
    'pier_corr08__merged_rep1_rep2.csv'
    ),
  na = ''
  )