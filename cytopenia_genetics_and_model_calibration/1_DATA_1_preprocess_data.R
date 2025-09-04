# Predict AML, MDS and MF

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Load library
source("mounts/research/src/Rfunctions/parallel_cluster.R")


######## Prepare data ########


## Read normalized healthy data
healthy = readRDS(paste0(import, "/lab_healthy_median.rds")) %>%
  dplyr::select(-n) %>%
  dplyr::rename(tulos_healthy = tulos)
## Healthy lab and demo
df_healthy_demo = readRDS(paste0(export, "/demo_healthy.rds")) %>%
  dplyr::distinct(henkilotunnus, sukupuoli_selite, syntymaaika_pvm, kuolinaika_pvm, .keep_all = FALSE)


# Lab tests to include
if (exists("labtests2")) { rm(labtests2) }
for (m in c("MDS", "MF", "de novo AML")) {  #"AML", "secondary AML"
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
  
  if (m == "MDS") {
    labtests2 = labtests1
  } else {
    labtests2 = rbind(labtests2, labtests1) %>%
      distinct(Variable, .keep_all = TRUE)
  }
}

# Loop
for (i in c("MDS", "MF", "de novo AML")) {
  
  print(paste0("Processing ", i))
  
  # Read demo
  i="MDS"
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  
  # Join data from healthy subjects
  df2_2 = df1 %>%
    dplyr::distinct()
  df2 = df2_2 %>%
    dplyr::left_join(healthy) %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy); rm(df2_2)
  
  # Process data
  df2 = df2 %>%
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
  
  
  lagged_data <- df2
  
  
  # Next create dynamic variables
  # Log
  writeLines(c(""), "mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt")
  sink("mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt", append=TRUE)
  
  # Export
  cat(paste("Export", length(unique(lagged_data$henkilotunnus)), "\n"))
  
  # Convert once to data.table
  setDT(lagged_data)
  
  # Pre-split by ID to avoid filtering inside the loop
  split_data <- split(lagged_data, by = "henkilotunnus")
  
  # Parallel backend
  
  # Counter for tracking
  lagged_data_1 = foreach(
    person_data = split_data,
    .packages = "data.table",
    .combine = rbind
  ) %dopar% {
    person_id <- unique(person_data$henkilotunnus)
    
    person_data <- person_data[order(time_to_dg)]
    
    person_data[, {
      res <- lapply(time_to_dg, function(x) {
        val_365  <- tulos_norm[time_to_dg < (x - 15) & time_to_dg < x & time_to_dg >= (x - 365)]
        val_1095 <- tulos_norm[time_to_dg < (x - 365) & time_to_dg >= (x - 1095)]
        val_1825 <- tulos_norm[time_to_dg < (x - 1095) & time_to_dg >= (x - 1825)]
        
        val_365_n = length(val_365)
        val_1095_n = length(val_1095)
        val_1825_n = length(val_1825)
        
        mean_365  <- if (length(val_365) > 0)  median(val_365, na.rm = TRUE) else NA_real_
        mean_1095 <- if (length(val_1095) > 0) median(val_1095, na.rm = TRUE) else NA_real_
        mean_1825 <- if (length(val_1825) > 0) median(val_1825, na.rm = TRUE) else NA_real_
        
        list(val_365_n = val_365_n, val_1095_n = val_1095_n, val_1825_n = val_1825_n,
             mean_365 = mean_365, mean_1095 = mean_1095, mean_1825 = mean_1825)
      })
      
      data.table(
        henkilotunnus = henkilotunnus,
        time_to_dg = time_to_dg,
        tulos_norm = tulos_norm,
        val_365_n = sapply(res, `[[`, "val_365_n"),
        val_1095_n = sapply(res, `[[`, "val_1095_n"),
        val_1825_n = sapply(res, `[[`, "val_1825_n"),
        mean_last_values_365_days = sapply(res, `[[`, "mean_365"),
        mean_last_values_1095_days = sapply(res, `[[`, "mean_1095"),
        mean_last_values_1825_days = sapply(res, `[[`, "mean_1825"),
        trend_from_365_d = tulos_norm - sapply(res, `[[`, "mean_365"),
        trend_from_1095_d = tulos_norm - sapply(res, `[[`, "mean_1095"),
        trend_from_1825_d = tulos_norm - sapply(res, `[[`, "mean_1825")
      )
    }, by = tutkimus_lyhenne_yksikko]
    
    
  }
  
  lagged_data_1 = as.data.frame(lagged_data_1) %>%
    dplyr::left_join(df2) %>%
    dplyr::arrange(henkilotunnus, naytteenotto_hetki, tutkimus_lyhenne) %>%
    dplyr::select(names(df2), starts_with("val"), starts_with("mean"), starts_with("trend"))
  
  sink(); sink()
  
  # Save
  lagged_data = lagged_data_1; rm(lagged_data_1)
  saveRDS(lagged_data, paste0(export, "/lagged_data_", i, ".rds"))
  fwrite(lagged_data, paste0(export, "/lagged_data_", i, ".csv"))
  
  print("Saving 2")
  
  
  ######## Reshape data ########
  
  
  lagged_data1 = lagged_data %>%
    dplyr::select(-c(naytteenotto_hetki, tulos, yksikko, viitearvot, tulos_viitearvo_ulkopuolella,
                     paatutkimus_lyhenne, syntymaaika_pvm, kuolinaika_pvm, dg_date_combined,
                     disease, time_to_dg_mo, tulos_healthy, tutkimus_lyhenne, age_group)) %>%
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
    dplyr::mutate(disease = i) %>%
    dplyr::filter(!(disease=="Healthy" & time_to_dg>-365)) %>%
    dplyr::mutate(event_1y = ifelse(time_to_dg >-365, -1/time_to_dg, 0))
  
  lagged_data1 = lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, everything())
  
  # Rename
  lagged_data1 = lagged_data1 %>%
    janitor::clean_names()
  
  
  # Remove variables
  names(lagged_data1)
  
  # Identify columns that contain both conditions
  cols_to_remove <- names(lagged_data1)[grepl("(ferrit|crp|trigly|kol|gf_re|p_tt|crea|b_la)", names(lagged_data1)) &
                                          grepl("(_d)|(tulos_norm)", names(lagged_data1))]
  # Remove those columns
  lagged_data1 <- lagged_data1[ , !(names(lagged_data1) %in% cols_to_remove)]
  
  
# Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "2.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "2.csv"))
  
  print("Saving 4")
  
  # Define predictor variables and outcome variable
  lagged_data1$disease = ifelse(lagged_data1$disease==i, 1, 0)
  unique(lagged_data1$disease)
  
  # Remove redundant variables
  lagged_data1 = lagged_data1 %>%
    dplyr::mutate(time_to_dg = -time_to_dg)
  
  # Calculate how many times lab tests have been taken within the last month
  lagged_data1 <- lagged_data1 %>%
    arrange(henkilotunnus, desc(time_to_dg)) %>%
    group_by(henkilotunnus) %>%
    mutate(rows_in_last_month = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 31 & time_to_dg >= t)),
           rows_in_last_year = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 366 & time_to_dg >= t))) %>%
    ungroup()
  
  # Remove data 5-10 years before diagnosis
  lagged_data1 <- lagged_data1 %>%
    dplyr::filter(disease == 0 | time_to_dg < 1824)
  
  # Remove atyply and i_ind if still there
  lagged_data1 = lagged_data1 %>%
    dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "i_ind", ignore.case = TRUE, names(lagged_data1))])) %>%
    dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "l_atyp_ly_percent", ignore.case = TRUE, names(lagged_data1))]))
  
  # Manually impute NA to 0.01 for erblast as these are by default 0 if not reported
  lagged_data1$b_erblast_e9_l_tulos_norm = ifelse(is.na(lagged_data1$b_erblast_e9_l_tulos_norm), 0.01, lagged_data1$b_erblast_e9_l_tulos_norm)
  lagged_data1$b_erblast_e9_l_trend_from_365_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_365_d), 0, lagged_data1$b_erblast_e9_l_trend_from_365_d)
  lagged_data1$b_erblast_e9_l_trend_from_1095_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1095_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1095_d)
  lagged_data1$b_erblast_e9_l_trend_from_1825_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1825_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1825_d)
  nrow(lagged_data1[!is.na(lagged_data1$b_monos_e9_l_tulos_norm),])
  
  # Exclude patients with wrong diagnosis
  exclude = readxl::read_xlsx(paste0("mounts/research/husdatalake/disease/general/exclude_pts.xlsx"))
  if (i == "de novo AML") {
    j = "AML"
  } else {
    j = i
  }
  lagged_data1 = lagged_data1 %>%
    dplyr::filter(!henkilotunnus %in% exclude[exclude$disease==j,]$henkilotunnus)
  
  # Remove secondary MF
  if (i == "MF") {
    pmf = readRDS("mounts/research/husdatalake/disease/processed_data/MF/secondary_MF.rds") %>%
      dplyr::filter(!disease == "Primary")
    lagged_data1 = lagged_data1 %>%
      dplyr::filter(!(disease == 1 & henkilotunnus %in% pmf$henkilotunnus))
  }
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "3.rds"))
  print("Saving 5")
  gc()
}

# Impute missing data
for (i in c("MDS", "MF", "de novo AML")) {
  
  print(i)
  
  # Load data
  lagged_data1 = readRDS(paste0(export, "/lagged_data_", i, "3.rds"))
  # Impute non NA
  ## Fill NAs with last non-NA
  ### First fill with linear interpolation
  lagged_data1 <- lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, rows_in_last_month, rows_in_last_year, contains("val_"), everything())
  lagged_data1 <- lagged_data1 %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.approx(., na.rm = FALSE)))
  ### Then fill upwards
  lagged_data1 <- lagged_data1 %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.locf(., fromLast = TRUE, na.rm = FALSE)))
  ### Then fill downwards
  lagged_data1 <- lagged_data1 %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.locf(., na.rm = FALSE)))
  ## For the rest (= all values missing for a lab test for a patient), impute 0
  lagged_data1 <- lagged_data1 %>%
    ungroup() %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):ncol(.)-1, ~ifelse(is.na(.), 0, .)))
  
  # Keep only adults
  lagged_data1 = lagged_data1 %>%
    dplyr::filter(age >= 18)

  
  # Log
  writeLines(c(""), "mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt")
  sink("mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt", append=TRUE)
  
  # Export
  cat(paste("Export", length(unique(lagged_data1$henkilotunnus)), "\n"))
  
  # Convert once to data.table
  setDT(lagged_data1)
  
  # Identify relevant columns
  norm_cols <- grep("_norm", names(lagged_data1), value = TRUE)
  
  # Pre-split by ID to avoid filtering inside the loop
  split_data <- split(lagged_data1, by = "henkilotunnus")
  
  # Parallelized loop
  lagged_data_1 <- foreach(
    person_data = split_data,
    .packages = "data.table",
    .combine = rbind
  ) %dopar% {
    
    # Ensure time is ordered
    person_data <- person_data[order(time_to_dg)]
    
    res_list <- vector("list", nrow(person_data))
    
    for (j in seq_len(nrow(person_data))) {
      x <- person_data$time_to_dg[j]
      
      out <- list(
        henkilotunnus = person_data$henkilotunnus[j],
        time_to_dg = x
      )
      
      for (col in norm_cols) {
        vals_365  <- person_data[[col]][person_data$time_to_dg < x & person_data$time_to_dg >= (x - 365)]
        vals_1095 <- person_data[[col]][person_data$time_to_dg < (x - 365) & person_data$time_to_dg >= (x - 1095)]
        vals_1825 <- person_data[[col]][person_data$time_to_dg < (x - 1095) & person_data$time_to_dg >= (x - 1825)]
        
        out[[paste0(col, "_val_365_n")]]  <- sum(!is.na(vals_365))
        out[[paste0(col, "_val_1095_n")]] <- sum(!is.na(vals_1095))
        out[[paste0(col, "_val_1825_n")]] <- sum(!is.na(vals_1825))
      }
      
      res_list[[j]] <- out
    }
    
    rbindlist(res_list)
  }
  
  lagged_data_1 = as.data.frame(lagged_data_1)
  
  lagged_data1 = lagged_data1 %>%
    dplyr::select(-contains("_val_")) %>%
    dplyr::left_join(lagged_data_1)
  lagged_data1 = lagged_data1 %>%
    dplyr::arrange(henkilotunnus, desc(time_to_dg)) %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y,
                  sukupuoli_selite, age, rows_in_last_month, rows_in_last_year,
                  sort(colnames(.)))
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/lagged_data_", i, "3.rds"))
  fwrite(lagged_data1, paste0(export, "/lagged_data_", i, "3.csv"))
}


##### HEALTHY #####


# Loop
filelist = list.files(path = "mounts/research/husdatalake/disease/data/Preleukemia/multilab3", pattern = "file[[:digit:]]*\\.parquet", full.names = TRUE)


for (i in 1:length(filelist)) {
  
  # for (i in c("secondary AML")) {
  print(paste0("Processing ", i))
  
  # Read demo
  df1 = arrow::read_parquet(filelist[i])
  
  # Read lab data for  healthy subjects
  df2_2 = df1 %>%
    dplyr::distinct()
  
  df2 = df2_2 %>%
    dplyr::left_join(healthy) %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy); rm(df2_2)
  
  # Process data
  df2 = df2 %>%
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
    # dplyr::filter(!henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::select(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age) %>% 
    distinct()
  dir.create(paste0(paste0(export, "/healthy")))
  fwrite(df3, paste0(export, "/healthy/lab_data_for_modelling_", i, ".csv"))
  
  print("Saving 1")
  
  lagged_data <- df2
  
  
  # Create dynamic variables
  # Log
  writeLines(c(""), "mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt")
  sink("mounts/research/husdatalake/disease/processed_data/Preleukemia/log_lagged_data.txt", append=TRUE)
  
  # Export
  cat(paste("Export", length(unique(lagged_data$henkilotunnus)), "\n"))
  
  # Convert once to data.table
  setDT(lagged_data)
  
  # Pre-split by ID to avoid filtering inside the loop
  split_data <- split(lagged_data, by = "henkilotunnus")
  
  # Counter for tracking
  lagged_data_1 = foreach(
    person_data = split_data,
    .packages = "data.table",
    .combine = rbind
  ) %dopar% {
    person_id <- unique(person_data$henkilotunnus)
    
    person_data <- person_data[order(time_to_dg)]
    
    person_data[, {
      res <- lapply(time_to_dg, function(x) {
        val_365  <- tulos_norm[time_to_dg < (x - 15) & time_to_dg < x & time_to_dg >= (x - 365)]
        val_1095 <- tulos_norm[time_to_dg < (x - 365) & time_to_dg >= (x - 1095)]
        val_1825 <- tulos_norm[time_to_dg < (x - 1095) & time_to_dg >= (x - 1825)]
        
        val_365_n = length(val_365)
        val_1095_n = length(val_1095)
        val_1825_n = length(val_1825)
        
        mean_365  <- if (length(val_365) > 0)  median(val_365, na.rm = TRUE) else NA_real_
        mean_1095 <- if (length(val_1095) > 0) median(val_1095, na.rm = TRUE) else NA_real_
        mean_1825 <- if (length(val_1825) > 0) median(val_1825, na.rm = TRUE) else NA_real_
        
        list(val_365_n = val_365_n, val_1095_n = val_1095_n, val_1825_n = val_1825_n,
             mean_365 = mean_365, mean_1095 = mean_1095, mean_1825 = mean_1825)
      })
      
      data.table(
        henkilotunnus = henkilotunnus,
        time_to_dg = time_to_dg,
        tulos_norm = tulos_norm,
        val_365_n = sapply(res, `[[`, "val_365_n"),
        val_1095_n = sapply(res, `[[`, "val_1095_n"),
        val_1825_n = sapply(res, `[[`, "val_1825_n"),
        mean_last_values_365_days = sapply(res, `[[`, "mean_365"),
        mean_last_values_1095_days = sapply(res, `[[`, "mean_1095"),
        mean_last_values_1825_days = sapply(res, `[[`, "mean_1825"),
        trend_from_365_d = tulos_norm - sapply(res, `[[`, "mean_365"),
        trend_from_1095_d = tulos_norm - sapply(res, `[[`, "mean_1095"),
        trend_from_1825_d = tulos_norm - sapply(res, `[[`, "mean_1825")
      )
    }, by = tutkimus_lyhenne_yksikko]
    
    
  }
  
  lagged_data_1 = as.data.frame(lagged_data_1) %>%
    dplyr::left_join(df2) %>%
    dplyr::arrange(henkilotunnus, naytteenotto_hetki, tutkimus_lyhenne) %>%
    dplyr::select(names(df2), starts_with("val"), starts_with("mean"), starts_with("trend"))
  
  sink(); sink()
  
  # Save
  lagged_data = lagged_data_1; rm(lagged_data_1)
  saveRDS(lagged_data, paste0(export, "/healthy/lagged_data_", i, "_.rds"))
  print("Saving 2")
  
  
  ######## Reshape data ########
  
  
  lagged_data1 = lagged_data %>%
    dplyr::select(-c(naytteenotto_hetki, tulos, yksikko, viitearvot, tulos_viitearvo_ulkopuolella,
                     paatutkimus_lyhenne, syntymaaika_pvm, kuolinaika_pvm, dg_date_combined,
                     disease, time_to_dg_mo, tulos_healthy, tutkimus_lyhenne, age_group)) %>%
    reshape2::melt(id.vars = c("henkilotunnus", "time_to_dg", "tutkimus_lyhenne_yksikko", "sukupuoli_selite", "age")) %>%
    dplyr::mutate(tutkimus_lyhenne_yksikko_var = paste0(tutkimus_lyhenne_yksikko, "_", variable)) %>%
    dplyr::select(-c(tutkimus_lyhenne_yksikko, variable)) %>%
    distinct()
  
  rm(df1, df2, lagged_data)
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
  
  # Identify columns that contain both conditions
  cols_to_remove <- names(lagged_data1)[grepl("(ferrit|crp|trigly|kol|gf_re|p_tt|crea|b_la)", names(lagged_data1)) &
                                          grepl("(_d)|(tulos_norm)", names(lagged_data1))]
  # Remove those columns
  lagged_data1 <- lagged_data1[ , !(names(lagged_data1) %in% cols_to_remove)]
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/healthy/lagged_data_", i, "_1.rds"))
  print("Saving 3")
  
  ## Process outcome for healthy
  lagged_data1 = lagged_data1 %>%
    dplyr::mutate(disease = "Healthy") %>%
    dplyr::filter(!(disease=="Healthy" & time_to_dg>-365)) %>%
    dplyr::mutate(event_1y = ifelse(time_to_dg >-365, -1/time_to_dg, 0))
  
  lagged_data1 = lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, everything())
  
  # Rename
  lagged_data1 = lagged_data1 %>%
    janitor::clean_names()
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/healthy/lagged_data_", i, "_2.rds"))
  
  print("Saving 4")
  
  # Define predictor variables and outcome variable
  lagged_data1$disease = 0
  
  # Remove redundant variables
  lagged_data1 = lagged_data1 %>%
    dplyr::mutate(time_to_dg = -time_to_dg)
  
  # Calculate how many times lab tests have been taken within the last month
  lagged_data1 <- lagged_data1 %>%
    arrange(henkilotunnus, desc(time_to_dg)) %>%
    group_by(henkilotunnus) %>%
    mutate(rows_in_last_month = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 31 & time_to_dg >= t)),
           rows_in_last_year = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 366 & time_to_dg >= t))) %>%
    ungroup()
  
  # Remove data 5-10 years before diagnosis
  lagged_data1 <- lagged_data1 %>%
    dplyr::filter(disease == 0 | time_to_dg < 1824)
  
  
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
  nrow(lagged_data1[!is.na(lagged_data1$b_monos_e9_l_tulos_norm),])
  
  
  # Remove if fewer than 3 rows per patient
  lagged_data1_sr = lagged_data1 %>%
    group_by(henkilotunnus) %>%
    summarise(n = n()) %>%
    dplyr::filter(n >= 2)
  lagged_data1 = lagged_data1 %>%
    dplyr::filter(henkilotunnus %in% lagged_data1_sr$henkilotunnus)
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/healthy/lagged_data_", i, "_3.rds"))
  print("Saving 5")
  
  gc()
}

# Impute missing data
listfiles1 = list.files(paste0(export, "/healthy"), pattern = "_3.rds", full.names = TRUE)
for (i in 1:length(listfiles1)) {
  print(i)
  # Load data
  lagged_data1 = readRDS(listfiles1[i])
  # Impute non NA
  ## Fill NAs with last non-NA
  ### First fill downwards with linear interpolation
  lagged_data1 <- lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, rows_in_last_month, rows_in_last_year, contains("val_"), everything())
  lagged_data1 <- lagged_data1 %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.approx(., na.rm = FALSE)))
  ### Then fill upwards
  lagged_data1 <- lagged_data1 %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.locf(., fromLast = TRUE, na.rm = FALSE)))
  ### Then fill downwards
  lagged_data1 <- lagged_data1 %>%
    group_by(henkilotunnus) %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):(ncol(.)-1), ~zoo::na.locf(., na.rm = FALSE)))
  ## For the rest (= all values missing for a lab test for a patient), impute 0
  lagged_data1 <- lagged_data1 %>%
    ungroup() %>%
    mutate(across((which(names(.)=="b_baso_e9_l_mean_last_values_1095_days")-1):ncol(.)-1, ~ifelse(is.na(.), 0, .)))
  
  # Keep only adults
  lagged_data1 = lagged_data1 %>%
    dplyr::filter(age >= 18)
  
  
  # Export
  print(paste("Export", length(unique(lagged_data1$henkilotunnus)), "\n"))
  
  # Convert once to data.table
  setDT(lagged_data1)
  
  # Identify relevant columns
  norm_cols <- grep("_norm", names(lagged_data1), value = TRUE)
  
  # Pre-split by ID to avoid filtering inside the loop
  split_data <- split(lagged_data1, by = "henkilotunnus")
  
  # Parallelized loop
  lagged_data_1 <- foreach(
    person_data = split_data,
    .packages = "data.table",
    .combine = rbind
  ) %dopar% {
    
    # Ensure time is ordered
    person_data <- person_data[order(time_to_dg)]
    
    res_list <- vector("list", nrow(person_data))
    
    for (j in seq_len(nrow(person_data))) {
      x <- person_data$time_to_dg[j]
      
      out <- list(
        henkilotunnus = person_data$henkilotunnus[j],
        time_to_dg = x
      )
      
      for (col in norm_cols) {
        vals_365  <- person_data[[col]][person_data$time_to_dg < x & person_data$time_to_dg >= (x - 365)]
        vals_1095 <- person_data[[col]][person_data$time_to_dg < (x - 365) & person_data$time_to_dg >= (x - 1095)]
        vals_1825 <- person_data[[col]][person_data$time_to_dg < (x - 1095) & person_data$time_to_dg >= (x - 1825)]
        
        out[[paste0(col, "_val_365_n")]]  <- sum(!is.na(vals_365))
        out[[paste0(col, "_val_1095_n")]] <- sum(!is.na(vals_1095))
        out[[paste0(col, "_val_1825_n")]] <- sum(!is.na(vals_1825))
      }
      
      res_list[[j]] <- out
    }
    
    rbindlist(res_list)
  }
  
  lagged_data_1 = as.data.frame(lagged_data_1)
  
  lagged_data1 = lagged_data1 %>%
    dplyr::select(-contains("_val_")) %>%
    dplyr::left_join(lagged_data_1)
  lagged_data1 = lagged_data1 %>%
    dplyr::arrange(henkilotunnus, desc(time_to_dg)) %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y,
                  sukupuoli_selite, age, rows_in_last_month, rows_in_last_year,
                  sort(colnames(.)))
  
  # Save
  saveRDS(lagged_data1, paste0(export, "/healthy/lagged_data_", i, "_4.rds"))
}

listfiles1 = list.files(paste0(export, "/healthy"), pattern = "_4.rds", full.names = TRUE)
lagged_data1 = data.frame()
for (i in 1:length(listfiles1)) {
  print(i)
  
  # Load data
  tmp = readRDS(listfiles1[i]) %>%
    # dplyr::select(-contains("_val_")) %>%
    dplyr::select(-contains("_mean_last_")) %>%
    dplyr::select(-ends_with("d_2")) %>%
    dplyr::select(-ends_with("norm_2")) %>%
    dplyr::select(-ends_with("_n_2"))
  
  # Identify columns that contain both conditions
  cols_to_remove <- names(tmp)[grepl("(ferrit|crp|trigly|kol|gf_re|p_tt|crea|b_la)", names(tmp)) &
                                          grepl("(_d)|(tulos_norm)", names(tmp))]
  # Remove those columns
  tmp <- tmp %>%
    dplyr::select(-any_of(cols_to_remove))
  
  lagged_data1 = bind_rows(lagged_data1, tmp)
}

# Export
saveRDS(lagged_data1, paste0(export, "/healthy/lagged_data_final.rds"))


##### Add disease data #####


# Read data

lagged_data1 = readRDS(paste0(export, "/healthy/lagged_data_final.rds"))
mf = readRDS(paste0(export, "/lagged_data_MF3.rds"))
saveRDS(mf, paste0(export, "/lagged_data_MF3_allpts.rds"))
mds = readRDS(paste0(export, "/lagged_data_MDS3.rds"))
saveRDS(mds, paste0(export, "/lagged_data_MDS3_allpts.rds"))
aml = readRDS(paste0(export, "/lagged_data_de novo AML3.rds"))
saveRDS(aml, paste0(export, "/lagged_data_de novo AML3_allpts.rds"))

# Remove if fewer than 3 rows per patient
lagged_data1_sr = mf %>%
  dplyr::filter(time_to_dg >= -10) %>%
  group_by(henkilotunnus) %>%
  summarise(n = n()) %>%
  dplyr::filter(n >= 2)
mf = mf %>%
  dplyr::filter(henkilotunnus %in% lagged_data1_sr$henkilotunnus)
lagged_data1_sr = mds %>%
  dplyr::filter(time_to_dg >= -10) %>%
  group_by(henkilotunnus) %>%
  summarise(n = n()) %>%
  dplyr::filter(n >= 2)
mds = mds %>%
  dplyr::filter(henkilotunnus %in% lagged_data1_sr$henkilotunnus)
lagged_data1_sr = aml %>%
  dplyr::filter(time_to_dg >= -10) %>%
  group_by(henkilotunnus) %>%
  summarise(n = n()) %>%
  dplyr::filter(n >= 2)
aml = aml %>%
  dplyr::filter(henkilotunnus %in% lagged_data1_sr$henkilotunnus)

# Keep only shared colnames
names1 = names(lagged_data1)
names_mf = names(mf)
names_mds = names(mds)
names_aml = names(aml)
shared_values <- Reduce(intersect, list(names1, names_mf, names_mds, names_aml))

# Combine
lagged_data1 = lagged_data1 %>%
  dplyr::select(all_of(shared_values))

# MF
mf = mf %>%
  dplyr::select(all_of(shared_values)) %>%
  bind_rows(lagged_data1)
fwrite(mf, paste0(export, "/lagged_data_MF3.csv"))
rm(mf)

# MDS
mds = mds %>%
  dplyr::select(all_of(shared_values)) %>%
  bind_rows(lagged_data1)
fwrite(mds, paste0(export, "/lagged_data_MDS3.csv"))
rm(mds)

# AML
aml = aml %>%
  dplyr::select(all_of(shared_values)) %>%
  bind_rows(lagged_data1)
fwrite(aml, paste0(export, "/lagged_data_de novo AML3.csv"))

# Edit lab_demo
## AML
df1 = readRDS(paste0(export, "/de novo AML_lab_demo.rds"))
aml = fread(paste0(export, "/lagged_data_de novo AML3.csv"), select = c("henkilotunnus", "disease")) %>% distinct()
aml = aml %>% dplyr::filter(disease == 1)
df1 = df1 %>%
  dplyr::filter(henkilotunnus %in% aml$henkilotunnus)
saveRDS(df1, paste0(export, "/de novo AML_lab_demo.rds"))
## MDS
df1 = readRDS(paste0(export, "/MDS_lab_demo.rds"))
aml = fread(paste0(export, "/lagged_data_MDS3.csv"), select = c("henkilotunnus", "disease")) %>% distinct()
aml = aml %>% dplyr::filter(disease == 1)
df1 = df1 %>%
  dplyr::filter(henkilotunnus %in% aml$henkilotunnus)
saveRDS(df1, paste0(export, "/MDS_lab_demo.rds"))
## MF
df1 = readRDS(paste0(export, "/MF_lab_demo.rds"))
aml = fread(paste0(export, "/lagged_data_MF3.csv"), select = c("henkilotunnus", "disease")) %>% distinct()
aml = aml %>% dplyr::filter(disease == 1)
df1 = df1 %>%
  dplyr::filter(henkilotunnus %in% aml$henkilotunnus)
saveRDS(df1, paste0(export, "/MF_lab_demo.rds"))
