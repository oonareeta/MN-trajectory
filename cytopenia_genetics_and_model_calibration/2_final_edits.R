# Last edits

# Author: Oscar Brück

# Libraries
source("./library.R")
library(doParallel)

# Create analysis cluster
nodes <- detectCores()-2
print(nodes)
cl <- makeCluster(nodes)
doParallel::registerDoParallel(nodes)

# General parameters
source("./parameters")


######## Prepare data ########


## Read normalized healthy data
healthy = readRDS(paste0(import, "/lab_healthy_median.rds")) %>%
  dplyr::select(-n) %>%
  dplyr::rename(tulos_healthy = tulos)
## Healthy lab and demo
df_healthy = readRDS(paste0(export, "/healthy_lab_demo.rds"))
df_healthy_demo = df_healthy %>%
  dplyr::distinct(henkilotunnus, sukupuoli_selite, syntymaaika_pvm, kuolinaika_pvm, .keep_all = FALSE)


# Lab tests to include
if (exists("labtests2")) { rm(labtests2) }
for (m in c("AML", "MDS", "MF", "de novo AML", "secondary AML")) {
  ## Adjust p
  ### Read linear regression tables - all
  labtests = readxl::read_xlsx(paste0(results, "/", m, "/linear_regression_results_5y.xlsx"))
  labtests$padj = p.adjust(labtests$P_Value, method = "BH")
  ### Read linear regression tables - 2 y pre-disease
  labtests1 = readxl::read_xlsx(paste0(results, "/", m, "/linear_regression_results_2y.xlsx"))
  labtests1$padj = p.adjust(labtests1$P_Value, method = "BH")
  ## Distinct tests
  labtests_0 = labtests %>%
    dplyr::filter(padj<0.05) %>%
    distinct(Variable, .keep_all = TRUE)
  labtests_1 = labtests1 %>%
    dplyr::filter(padj<0.05) %>%
    distinct(Variable, .keep_all = TRUE)
  colnames(labtests_1)[2:4] = paste0(colnames(labtests_1)[2:4], "_2y")
  ## Keep only variables predictive both 5-0 and 2-0 years prior to disease
  labtests1 = inner_join(labtests_0, labtests_1)
  
  if (m == "AML") {
    labtests2 = labtests1
  } else {
    labtests2 = rbind(labtests2, labtests1) %>%
      distinct(Variable, .keep_all = TRUE)
  }
}
# Remove unnecessary variables
labtests = labtests %>%
  dplyr::filter(!str_detect(Variable, "(?i)(GHb-A1C|I-ind|-Kol|Trigly|Gluc|HbA1c|atyply|Lacta|Folaa|Trfesa|Transf|fP-Fe)"))


# Loop
for (i in c("MDS", "MF", "de novo AML")) {
  
  # for (i in c("secondary AML")) {
  print(paste0("Processing ", i))
  
  # Read demo
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read lab data for 10000 healthy subjects
  df = readRDS(paste0(export, "/healthy_lab_demo_all_10000subjects_post_", i, ".rds")) ## Post-processed
  df3 = readRDS(paste0(export, "/healthy_lab_demo_all_10000subjects.rds")) %>% ## Pre-processed with all variables
    dplyr::distinct(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tutkimus_lyhenne, age, age_group, sukupuoli_selite, .keep_all=FALSE)
  
  # Impact of age and gender
  df2_1 = df %>%
    dplyr::inner_join(df3, relationship = "many-to-many")
  df2_2 = df %>%
    dplyr::inner_join(df1 %>%
                        dplyr::select(names(df3)) %>%
                        dplyr::distinct(), relationship = "many-to-many")
  df2 = rbind(df2_1, df2_2) %>%
    dplyr::left_join(healthy) %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy); rm(df2_1, df2_2)
  
  
  # Process data
  df2 = df2 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% labtests2$Variable) %>%
    arrange(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg)
  
  
  
  # If two values on the same day AND the test names are starting with "L-" AND the other is 0.01 > remove the 0.01
  tmp = df2 %>%
    ungroup() %>%
    dplyr::filter(str_detect(tutkimus_lyhenne_yksikko, "^L-")) %>%
    group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n > 1)
  tmp1 = df2 %>%
    dplyr::ungroup() %>%
    dplyr::filter(str_detect(tutkimus_lyhenne_yksikko, "^L-")) %>%
    dplyr::left_join(tmp) %>%
    dplyr::select(-n) %>%
    dplyr::arrange(time_to_dg, desc(tulos)) %>%
    dplyr::group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg) %>%
    dplyr::slice(1)
  df2 = df2 %>%
    dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^L-")) %>%
    dplyr::full_join(tmp1)
  
  # Export data
  df3 = df2 %>%
    dplyr::select(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age) %>% 
    distinct()
  fwrite(df3, paste0(export, "/lab_data_for_modelling_", i, ".csv"))
  
  print("Saving 1")
  
  
  # As we have run already the data for healthy subjects, we will run only disease data for MF and de novo AML
  # In the end of the script we will add the healthy data
  if (i %in% c("MF", "de novo AML")) {
    df2 = df2 %>%
      dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus))
  }
  lagged_data <- df2
  
  
  # Log
  writeLines(c(""), "mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt")
  sink("mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt", append=TRUE)
  
  # Export
  cat(paste("Export", length(unique(lagged_data$henkilotunnus)), "\n"))
  
  lagged_data_1 = foreach(i = unique(lagged_data$henkilotunnus),.combine=bind_rows) %dopar% {
    cat(paste(i, "\n"))
    
    lagged_data %>%
      ungroup() %>%
      dplyr::filter(henkilotunnus == i) %>%
      group_by(tutkimus_lyhenne_yksikko) %>% #henkilotunnus,
      mutate(
        # Find the value from the last observation 1-1800 days before sampling
        last_values_365_days = map_dbl(time_to_dg, function(x) {
          match_row <- cur_data() %>%
            filter(time_to_dg < x & time_to_dg >= (x - 365)) %>%
            slice_min(time_to_dg, n = 1, with_ties = FALSE)
          if (nrow(match_row) > 0) match_row$tulos_norm else NA_real_
        }),
        # Calculate the median of the period ≥ 365 days ago
        mean_last_values_365_days = median(last_values_365_days, na.rm=TRUE),
        # Calculate the difference from the current value to the mean of the last three values ≥ 365 days ago
        trend_from_365_d = tulos_norm - mean_last_values_365_days,
        # Find the value from the last observation 365-1095 days before sampling
        last_values_1095_days = map_dbl(time_to_dg, function(x) {
          match_row <- cur_data() %>%
            filter(time_to_dg <= (x - 365) & time_to_dg >= (x - 1095)) %>%
            slice_min(time_to_dg, n = 1, with_ties = FALSE)
          if (nrow(match_row) > 0) match_row$tulos_norm else NA_real_
        }),
        # Calculate the median of the period 365-1095 days before sampling
        mean_last_values_1095_days = median(last_values_1095_days, na.rm=TRUE),
        # Calculate the difference from the current value to the mean of the last three values 365-1095 days before sampling
        trend_from_1095_d = tulos_norm - mean_last_values_1095_days,
        # Find the value from the last observation 365-1095 days before sampling
        last_values_1825_days = map_dbl(time_to_dg, function(x) {
          match_row <- cur_data() %>%
            filter(time_to_dg <= (x - 1095) & time_to_dg >= (x - 1825)) %>%
            slice_min(time_to_dg, n = 1, with_ties = FALSE)
          if (nrow(match_row) > 0) match_row$tulos_norm else NA_real_
        }),
        # Calculate the median of the period 1095-1825 days before sampling
        mean_last_values_1825_days = median(last_values_1825_days, na.rm=TRUE),
        # Calculate the difference from the current value to the mean of the last three values 1095-1825 days before sampling
        trend_from_1825_d = tulos_norm - mean_last_values_1825_days,
      ) %>%
    ungroup() %>%
      dplyr::select(-c(last_values_365_days, last_values_1095_days, last_values_1825_days))
  }
  sink(); sink()
  
  # Save
  lagged_data = lagged_data_1; rm(lagged_data_1)
  saveRDS(lagged_data, paste0(export, "/lagged_data_", i, ".rds"))
  fwrite(lagged_data, paste0(export, "/lagged_data_", i, ".csv"))
  print("Saving 2")
  
  
  ######## Reshape data ########
  
  
  lagged_data1 = lagged_data %>%
    dplyr::select(-c(tulos, tulos_healthy, tutkimus_lyhenne, age_group)) %>%
    ## We will need to remove blast proportion as this is indicative of leukemia
    reshape2::melt(id.vars = c("henkilotunnus", "time_to_dg", "tutkimus_lyhenne_yksikko", "sukupuoli_selite", "age")) %>%
    dplyr::mutate(tutkimus_lyhenne_yksikko_var = paste0(tutkimus_lyhenne_yksikko, "_", variable)) %>%
    dplyr::select(-c(tutkimus_lyhenne_yksikko, variable)) %>%
    distinct()
  
  rm(df, df1, df2, lagged_data)
  gc()
  
  # Export
  ptsn = length(unique(lagged_data1$henkilotunnus))
  ptsl = unique(lagged_data1$henkilotunnus)
  
  tmp = data.frame()
  for (jj in (1:ptsn)) {
    print(paste0(jj, "/", ptsn))
    
    tmp = bind_rows(tmp, lagged_data1 %>%
                      dplyr::filter(henkilotunnus == ptsl[jj]) %>%
                      reshape2::dcast(henkilotunnus+time_to_dg+sukupuoli_selite+age~tutkimus_lyhenne_yksikko_var,
                                      fun.aggregate = function(x) mean(x, na.rm = TRUE)))
    
  }
  lagged_data1 = tmp; rm(tmp)
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "1.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "1.csv"))
  
  print("Saving 3")
  
  ## Process outcome for healthy
  lagged_data1 = lagged_data1 %>%
    dplyr::left_join(df_healthy %>% dplyr::distinct(henkilotunnus) %>% dplyr::mutate(disease="Healthy")) %>%
    dplyr::mutate(disease = ifelse(is.na(disease), i, disease)) %>%
    dplyr::filter(!(disease=="Healthy" & time_to_dg>-365)) %>%
    dplyr::mutate(event_1y = ifelse(time_to_dg >-365, -1/time_to_dg, 0))
  lagged_data1 = lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, everything())
  # Rename
  lagged_data1 = lagged_data1 %>%
    janitor::clean_names()
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "2.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "2.csv"))
  # Add the remaining healthy as you run the healthy data in two batches (only with MDS)
  if (i == "MDS") {
    lagged_data1 = rbind(lagged_data1, readRDS(paste0(export, "/lagged_data_Healthy2.rds")))
  }
  
  print("Saving 4")
  
  # Define predictor variables and outcome variable
  predictors1 <- colnames(lagged_data1)[7:(ncol(lagged_data1))]
  lagged_data1$disease = ifelse(lagged_data1$disease==i, 1, 0)
  
  # Remove redundant variables
  lagged_data1 = lagged_data1 %>%
    dplyr::mutate(time_to_dg = -time_to_dg)
  
  # Calculate how many times lab tests have been taken within the last month
  lagged_data1 <- lagged_data1 %>%
    arrange(henkilotunnus, desc(time_to_dg)) %>%
    group_by(henkilotunnus) %>%
    mutate(rows_in_last_month = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 31 & time_to_dg >= t))) %>%
    ungroup()
  
  # Remove data 5-10 years before diagnosis
  lagged_data1 <- lagged_data1 %>%
    dplyr::filter(disease == 0 | time_to_dg < 1824)
  
  
  # Add healthy data
  if (i != "MDS") {
    tmp = readRDS(paste0(export, "/lagged_data_MDS3.rds")) %>%
      dplyr::filter(disease == 0)
    lagged_data1 = bind_rows(lagged_data1, tmp)
  }
  
  # Remove atyply and i_ind if still there
  lagged_data1 = lagged_data1 %>%
    dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "i_ind", ignore.case = TRUE, names(lagged_data1))])) %>%
    dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "l_atyp_ly_percent", ignore.case = TRUE, names(lagged_data1))]))
  
  # Manually impute NA to 0.01 for erblast as these are by default 0 if not reported
  lagged_data1$b_erblast_e9_l_tulos_norm = ifelse(is.na(lagged_data1$b_erblast_e9_l_tulos_norm), 0.01, lagged_data1$b_erblast_e9_l_tulos_norm)
  lagged_data1$b_erblast_e9_l_mean_last_values_365_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_365_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_365_days)
  lagged_data1$b_erblast_e9_l_mean_last_values_1825_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_1825_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_1825_days)
  lagged_data1$b_erblast_e9_l_mean_last_values_1095_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_1095_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_1095_days)
  lagged_data1$b_erblast_e9_l_trend_from_365_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_365_d), 0, lagged_data1$b_erblast_e9_l_trend_from_365_d)
  lagged_data1$b_erblast_e9_l_trend_from_1095_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1095_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1095_d)
  lagged_data1$b_erblast_e9_l_trend_from_1825_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1825_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1825_d)
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "3.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "3.csv"))
  print("Saving 5")
  gc()
}

# Impute missing data
# # alaspäin = potilaskohtainen lineaarinen interpolaatio (ei ollutkaan niin hidas, joku 2 min)
# # ylöspäin = potilaskohtainen viimeisin nonNA
# # jos kaikki on NA = koko aineiston mediaani

for (i in c("MDS", "MF", "de novo AML")) {
  
  # Load data
  lagged_data1 = readRDS(paste0(export, "/lagged_data_", i, "3.rds"))
  # Impute non NA
  ## Fill NAs with last non-NA
  ### First fill downwards with linear interpolation
  lagged_data1 <- lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, rows_in_last_month, everything())
  lagged_data1 <- lagged_data1 %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    group_by(henkilotunnus) %>%
    mutate(across(6:(ncol(.)-1), ~zoo::na.approx(., na.rm = FALSE)))
  ### Then fill upwards
  lagged_data1 <- lagged_data1 %>%
    group_by(henkilotunnus) %>%
    mutate(across(6:(ncol(.)-1), ~zoo::na.locf(., fromLast = TRUE, na.rm = FALSE)))
  ## For the rest, impute median by group
  lagged_data1 <- lagged_data1 %>%
    ungroup() %>%
    mutate(across(6:ncol(.)-1, ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "3.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "3.csv"))
}

