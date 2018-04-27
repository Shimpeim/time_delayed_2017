result_1 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_1_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_2 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_2_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_3 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_3_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_4 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_4_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_5 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_5_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_6 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_6_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_7 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_7_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_8 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_8_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_9 <- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_9_all_permuted_Rho_Q.val_observed.csv',sep=';')
result_10<- read.csv('./Simulation_Result/alpha_error_simulation_permtest_itt__permute_10000_batch_10_all_permuted_Rho_Q.val_observed.csv',sep=';')

result <- result_1 %>%
  rbind(result_2) %>%
  rbind(result_3) %>%
  rbind(result_4) %>%
  rbind(result_5) %>%
  rbind(result_6) %>%
  rbind(result_7) %>%
  rbind(result_8) %>%
  rbind(result_9) %>%
  rbind(result_10) %>%
  summarize(
    pval.pt = sum(Flg_pval.pt),
    pval.sp = sum(Flg_pval.sp)
  )