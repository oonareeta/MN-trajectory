# Process data

# Author: Oscar Brück

# Libraries
source("./library.R")
library(broom)

# General parameters
source("./parameters")


################ Load data #################################################################


# Read data
df = readRDS(paste0(export, "/lab_aml_mds_mf_corrected.rds"))
a = readRDS(paste0(export, "/CBC_tests_to_include.rds"))


# Exclude patients with wrong diagnosis
exclude = readxl::read_xlsx(paste0("mounts/research/husdatalake/disease/general/exclude_pts.xlsx"))
# if (i == "de novo AML") {
#   j = "AML"
# } else {
#   j = i
# }
df = df %>%
  dplyr::filter(!(disease == "AML" & henkilotunnus %in% exclude[exclude$disease=="AML",]$henkilotunnus))
df = df %>%
  dplyr::filter(!(disease == "MF" & henkilotunnus %in% exclude[exclude$disease=="MF",]$henkilotunnus))
df = df %>%
  dplyr::filter(!(disease == "MDS" & henkilotunnus %in% exclude[exclude$disease=="MDS",]$henkilotunnus))

# Remove secondary MF
pmf = readRDS("mounts/research/husdatalake/disease/processed_data/MF/secondary_MF.rds") %>%
  dplyr::filter(!disease == "Primary")
df = df %>%
  dplyr::filter(!(disease == "MF" & henkilotunnus %in% pmf$henkilotunnus))



# Time to dg
df = df %>%
  dplyr::mutate(
    time_to_dg = dg_date_combined - naytteenottoaika,
    time_to_dg = as.numeric(time_to_dg)
  ) %>%
  ## Filter tests without digits
  dplyr::filter(
    str_detect(tulos, "[[:digit:]]")
  ) %>%
  mutate(
    tulos = as.numeric(tulos)
  ) %>%
  dplyr::filter(time_to_dg < (10*365) & time_to_dg > 30)


# De novo AML
df_mds_mf = df %>%
  dplyr::filter(disease %in% c("MDS", "MF"))
df_aml_denovo = df %>%
  dplyr::filter(disease == "AML") %>%
  dplyr::filter(!henkilotunnus %in% df_mds_mf$henkilotunnus) %>%
  dplyr::mutate(disease = "de novo AML")
df_aml_sec = df %>%
  dplyr::filter(disease == "AML") %>%
  dplyr::filter(henkilotunnus %in% df_mds_mf$henkilotunnus) %>%
  dplyr::mutate(disease = "secondary AML")
df = df %>%
  bind_rows(df_aml_denovo) %>%
  bind_rows(df_aml_sec)


# Remove tests ordered rarely
## Summarise
df_n_aml1 = df %>%
  dplyr::filter(disease=="AML") %>%
  distinct(henkilotunnus, tutkimus_lyhenne) %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 100)
df_n_aml2 = df %>%
  dplyr::filter(disease=="AML") %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 300)
df_n_mds1 = df %>%
  dplyr::filter(disease=="MDS") %>%
  distinct(henkilotunnus, tutkimus_lyhenne) %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 100)
df_n_mds2 = df %>%
  dplyr::filter(disease=="MDS") %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 300)
df_n_mf1 = df %>%
  dplyr::filter(disease=="MF") %>%
  distinct(henkilotunnus, tutkimus_lyhenne) %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 50)
df_n_mf2 = df %>%
  dplyr::filter(disease=="MF") %>%
  group_by(tutkimus_lyhenne) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 150)
df_lab_tests = rbind(df_n_aml1, df_n_aml2, df_n_mds1, df_n_mds2, df_n_mf1, df_n_mf2) %>%
  distinct(tutkimus_lyhenne)
df_negative = df %>%
  dplyr::filter(tulos<0) %>%
  distinct(tutkimus_lyhenne)
saveRDS(df_lab_tests, paste0(export, "/df_lab_tests.rds"))
saveRDS(df_negative, paste0(export, "/df_negative.rds"))


# Loop over distinct diseases
for (i in unique(df$disease)) {
  print(i)
  df1 = df %>%
    dplyr::filter(tutkimus_lyhenne %in% a$Var1 |
                    (tutkimus_lyhenne %in% df_lab_tests$tutkimus_lyhenne)) %>%
    dplyr::filter(disease == i) %>%
    dplyr::filter(!tutkimus_lyhenne %in% df_negative$tutkimus_lyhenne) %>%
    dplyr::filter(!is.na(tulos)) %>%
    dplyr::mutate(tulos = as.numeric(tulos),
                  tulos = tulos + 0.01,
                  time_to_dg_mo = ceiling(time_to_dg/365.24),
                  time_to_dg = -time_to_dg,
                  time_to_dg_mo = -time_to_dg_mo,
                  age = round(((naytteenottoaika-as.Date(syntymaaika_pvm))/365.24), 1),
                  tutkimus_lyhenne_yksikko = gsub("Ion\\.", "Ion", gsub(" ", "", tutkimus_lyhenne)),
                  tutkimus_lyhenne_yksikko = ifelse(str_detect(tutkimus_lyhenne_yksikko, "^L") & yksikko == "", paste0(tutkimus_lyhenne_yksikko, " (%)"),
                                                    ifelse(yksikko == "", tutkimus_lyhenne_yksikko, paste0(tutkimus_lyhenne_yksikko, " (", gsub("_", " ", yksikko), ")"))),
                  tutkimus_lyhenne_yksikko = gsub("Liuskat", "Segment",
                                                  gsub("Retik", "Retic",
                                                       gsub("Lymf", "Ly",
                                                            gsub("Punas", "RBC",
                                                                 gsub("Sauvatn", "Band",
                                                                      gsub("-Kj", " Conj",
                                                                           gsub("Gluk", "Gluc",
                                                                                gsub("Krea", "Crea",
                                                                                     gsub("Lier", "Cast",
                                                                                          gsub("AtyypLymf", " AtypLy",
                                                                                               gsub("Laktaat", "Lactate",
                                                                                                    gsub("Leuk", "WBC",
                                                                                                         gsub("Trom", "PLT", tutkimus_lyhenne_yksikko)))))))))))))) %>%
    dplyr::filter(!tutkimus_lyhenne_yksikko  %in% c("P-PSA (ug/l)", "P-PSA-V (ug/l)", "P-TnI (ug/l)", "P-TnT (ug/l)", "S-ANA", "S-B12-Vit (pmol/l)"))
  
  # Define the bin breaks (including rightmost bin for ages 100 and above)
  breaks <- c(seq(0, 100, by = 5), 200)
  # Label the bins with descriptive text
  bin_labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])
  bin_labels = bin_labels[-21]
  bin_labels[length(bin_labels)+1] <- "100+"  # Add label for the last bin
  # Create a new column named "age_group" with binned categories
  df1$age_group <- cut(as.numeric(df1$age), breaks = breaks, right = FALSE, labels = bin_labels)
  
  ## Save data
  saveRDS(df1, paste0(export, "/", i, "_lab_demo.rds"))
  
  
  # Fit linear regression models for each variable
  tutkimus_lyhenne_yksikko_list = df1 %>%
    dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^(Pt|S-|P-|aB-|vB-|U-|-I-|fS|fP|fE|B-HbA1c|B-GHb|B-Reti)")) %>%
    dplyr::filter(!tutkimus_lyhenne_yksikko  %in% c("P-PSA (ug/l)", "P-PSA-V (ug/l)", "P-TnI (ug/l)", "P-TnT (ug/l)",
                                                    "S-ANA", "S-B12-Vit (pmol/l)")) %>%
    distinct(tutkimus_lyhenne_yksikko)
  
  
  ##### Impact of age and gender ##### 
  
  
  # Read data on healthy
  healthy = readRDS(paste0(import, "/lab_healthy_median.rds")) %>%
    dplyr::select(-n) %>%
    dplyr::rename(tulos_healthy = tulos)
  df1 = df1 %>%
    dplyr::left_join(healthy)
  df1 = df1 %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy)
  
  
  ##### LINEAR REGRESSION ##### 
  
  
  # Linear regression with only the last 5 years. Repeat for the last 4, 3, 2 and 1 year. Save the results for a later heatmap.
  
  
  # 5 years
  # Exclude variables with no data 5 years prior to the event
  df1_tmp = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) %>%
    dplyr::filter(time_to_dg >= (-5*365.24))
  vb = as.data.frame(table(df1_tmp$tutkimus_lyhenne_yksikko, df1_tmp$time_to_dg_mo)) %>%
    dplyr::mutate(Var2 = as.numeric(as.character(Var2)),
                  Freq = as.numeric(as.character(Freq))) %>%
    dplyr::filter(Var2 >= (-5) & Freq > 0) %>%  # Var2 = years before dg
    distinct(Var1)
  tt = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko[tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko %in% vb$Var1]
  
  # Linear regression
  models <- lapply(tt, function(variable) {
    lm(tulos_norm ~ time_to_dg, data = filter(df1[df1$time_to_dg >= (-5*365.24),], tutkimus_lyhenne_yksikko == variable)) # df1$time_to_dg <= (-6*(365.24/12)) & 
  })
  print(paste0("The number of patients with lab data before diagnosis is ", length(unique(df1[df1$time_to_dg_mo<0,]$henkilotunnus))))
  
  # Extract coefficients and p-values
  coefficients <- lapply(models, coef)
  # p_values <- lapply(models, function(model) summary(model)$coefficients[, 4])
  p_values <- lapply(models, function(model) summary(model)$coefficients[, 4])
  # Combine coefficients and p-values into a dataframe
  results1 <- data.frame(
    Variable = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko,
    Coefficient = sapply(coefficients, "[[", 2),
    P_Value = sapply(p_values, "[[", 2)
  ) %>%
    dplyr::arrange(P_Value)
  ## Save
  dir.create(paste0(results, "/", i, "/lab_boxplots"), recursive = TRUE)
  writexl::write_xlsx(results1, paste0(results, "/", i, "/linear_regression_results_5y.xlsx"))
  
  
  
  # 4 years
  # Exclude variables with no data 4 years prior to the event
  df1_tmp = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) %>%
    dplyr::filter(time_to_dg >= (-4*365.24))
  vb = as.data.frame(table(df1_tmp$tutkimus_lyhenne_yksikko, df1_tmp$time_to_dg_mo)) %>%
    dplyr::mutate(Var2 = as.numeric(as.character(Var2)),
                  Freq = as.numeric(as.character(Freq))) %>%
    dplyr::filter(Var2 >= (-4) & Freq > 0) %>%  # Var2 = 0.5 years before dg
    distinct(Var1)
  tt = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko[tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko %in% vb$Var1]
  
  models_4y <- lapply(tt, function(variable) {
    lm(tulos_norm ~ time_to_dg, data = filter(df1[df1$time_to_dg >= (-4*365.24),], tutkimus_lyhenne_yksikko == variable))  #df1$time_to_dg <= (-6*(365.24/12)) &
  })
  
  
  # Extract coefficients and p-values
  coefficients <- lapply(models_4y, coef)
  p_values <- lapply(models_4y, function(model) summary(model)$coefficients[, 4])
  
  # Combine coefficients and p-values into a dataframe
  results1_4y <- data.frame(
    Variable = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko,
    Coefficient = sapply(coefficients, "[[", 2),
    P_Value = sapply(p_values, "[[", 2)
  ) %>%
    dplyr::arrange(P_Value)
  ## Save
  writexl::write_xlsx(results1_4y, paste0(results, "/", i, "/linear_regression_results_4y.xlsx"))
  
  
  # 3 years
  # Exclude variables with no data 3 years prior to the event
  df1_tmp = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) %>%
    dplyr::filter(time_to_dg >= (-3*365.24))
  vb = as.data.frame(table(df1_tmp$tutkimus_lyhenne_yksikko, df1_tmp$time_to_dg_mo)) %>%
    dplyr::mutate(Var2 = as.numeric(as.character(Var2)),
                  Freq = as.numeric(as.character(Freq))) %>%
    dplyr::filter(Var2 >= (-3) & Freq > 0) %>%  # Var2 = 0.5 years before dg
    distinct(Var1)
  tt = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko[tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko %in% vb$Var1]
  
  models_3y <- lapply(tt, function(variable) {
    lm(tulos_norm ~ time_to_dg, data = filter(df1[df1$time_to_dg >= (-3*365.24),], tutkimus_lyhenne_yksikko == variable))  #df1$time_to_dg <= (-6*(365.24/12)) &
  })
  
  
  # Extract coefficients and p-values
  coefficients <- lapply(models_3y, coef)
  p_values <- lapply(models_3y, function(model) summary(model)$coefficients[, 4])
  
  # Combine coefficients and p-values into a dataframe
  results1_3y <- data.frame(
    Variable = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko,
    Coefficient = sapply(coefficients, "[[", 2),
    P_Value = sapply(p_values, "[[", 2)
  ) %>%
    dplyr::arrange(P_Value)
  ## Save
  writexl::write_xlsx(results1_3y, paste0(results, "/", i, "/linear_regression_results_3y.xlsx"))
  
  
  # 2 years
  # Exclude variables with no data 2 years prior to the event
  df1_tmp = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) %>%
    dplyr::filter(time_to_dg >= (-2*365.24))
  vb = as.data.frame(table(df1_tmp$tutkimus_lyhenne_yksikko, df1_tmp$time_to_dg_mo)) %>%
    dplyr::mutate(Var2 = as.numeric(as.character(Var2)),
                  Freq = as.numeric(as.character(Freq))) %>%
    dplyr::filter(Var2 >= (-2) & Freq > 0) %>%  # Var2 = 0.5 years before dg
    distinct(Var1)
  tt = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko[tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko %in% vb$Var1]
  
  models_2y <- lapply(tt, function(variable) {
    lm(tulos_norm ~ time_to_dg, data = filter(df1[df1$time_to_dg >= (-2*365.24),], tutkimus_lyhenne_yksikko == variable))  #df1$time_to_dg <= (-6*(365.24/12)) &
  })
  
  
  # Extract coefficients and p-values
  coefficients <- lapply(models_2y, coef)
  p_values <- lapply(models_2y, function(model) summary(model)$coefficients[, 4])
  
  # Combine coefficients and p-values into a dataframe
  results1_2y <- data.frame(
    Variable = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko,
    Coefficient = sapply(coefficients, "[[", 2),
    P_Value = sapply(p_values, "[[", 2)
  ) %>%
    dplyr::arrange(P_Value)
  ## Save
  writexl::write_xlsx(results1_2y, paste0(results, "/", i, "/linear_regression_results_2y.xlsx"))
  
  
  # 1 years
  # Exclude variables with no data 1 years prior to the event
  df1_tmp = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) %>%
    dplyr::filter(time_to_dg >= (-1*365.24))
  vb = as.data.frame(table(df1_tmp$tutkimus_lyhenne_yksikko, df1_tmp$time_to_dg_mo)) %>%
    dplyr::mutate(Var2 = as.numeric(as.character(Var2)),
                  Freq = as.numeric(as.character(Freq))) %>%
    dplyr::filter(Var2 >= (-1) & Freq > 0) %>%  # Var2 = 0.5 years before dg
    distinct(Var1)
  tt = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko[tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko %in% vb$Var1]
  
  models_1y <- lapply(tt, function(variable) {
    lm(tulos_norm ~ time_to_dg, data = filter(df1[df1$time_to_dg >= (-1*365.24),], tutkimus_lyhenne_yksikko == variable))  #df1$time_to_dg <= (-6*(365.24/12)) &
  })
  
  
  # Extract coefficients and p-values
  coefficients <- lapply(models_1y, coef)
  p_values <- lapply(models_1y, function(model) summary(model)$coefficients[, 4])
  
  # Combine coefficients and p-values into a dataframe
  results1_1y <- data.frame(
    Variable = tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko,
    Coefficient = sapply(coefficients, "[[", 2),
    P_Value = sapply(p_values, "[[", 2)
  ) %>%
    dplyr::arrange(P_Value)
  ## Save
  writexl::write_xlsx(results1_1y, paste0(results, "/", i, "/linear_regression_results_1y.xlsx"))
  
}


##### PLOT #####


for (i in c("MDS", "de novo AML", "MF")) {
  
  print(i)
  
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read data on healthy
  healthy = readRDS(paste0(import, "/lab_healthy_median.rds")) %>%
    dplyr::select(-n) %>%
    dplyr::rename(tulos_healthy = tulos)
  df1 = df1 %>%
    dplyr::left_join(healthy)
  df1 = df1 %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy)
  
  # La to ESR
  df1 = df1 %>%
    dplyr::mutate(tutkimus_lyhenne_yksikko = ifelse(tutkimus_lyhenne_yksikko=="B-La (mm/h)", "B-ESR (mm/h)", tutkimus_lyhenne_yksikko))
  
  # Distinct tests
  tutkimus_lyhenne_yksikko_list = df1 %>%
    distinct(tutkimus_lyhenne_yksikko)
  
  # Colors
  if (i == "MDS") {
    cols1 = "#348ABD"
  } else if (i == "de novo AML") {
    cols1 = "#BF9F45"
  } else if (i == "MF") {
    cols1 = "#2B6E2A"
  } 
  
  # Summarise median by patient and laboratory test
  df1 = df1 %>%
    dplyr::group_by(henkilotunnus, time_to_dg_mo, tutkimus_lyhenne_yksikko) %>%
    mutate(tulos_norm = mean(tulos_norm, na.rm=TRUE),
           tulos = mean(tulos, na.rm=TRUE)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Plot each regression with original values
  print("original")
  for (j in tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) {
    g = ggplot(df1[df1$time_to_dg_mo > -5.1 & df1$tutkimus_lyhenne_yksikko == j,], aes(x = time_to_dg_mo, y = tulos, group=time_to_dg_mo)) +
      geom_jitter(width = 0.2, size = 0.2, color = "black") +
      geom_boxplot(alpha = 0.3, fill = cols1, alpha = 0.3, color = "black", outliers = FALSE) +
      geom_hline(yintercept = median(df1[df1$tutkimus_lyhenne_yksikko == j & df1$time_to_dg_mo == min(df1[df1$time_to_dg_mo > -5.1,]$time_to_dg_mo, na.rm = TRUE),]$tulos, na.rm = TRUE),
                 size = 1, alpha = 1, color = "#e41a1c", lty="11") +
      labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
           y=gsub("^L-", "", j), x=paste0("Time to diagnosis (years)")) +
      scale_x_continuous(breaks = seq(from = -5, to = -1, by = 1), labels = seq(from = -5, to = -1, by = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title = element_text(size=14, colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size=14, face="bold", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none")
    g
    # Export
    ggsave(plot = g,
           filename = paste0(results, "/", i, "/lab_boxplots/", janitor::make_clean_names(j), ".png"),
           height = 4, width = 4, units = "in", dpi = 300)
  }
  
  # Plot each regression with normalized values
  print("norm")
  
  for (j in tutkimus_lyhenne_yksikko_list$tutkimus_lyhenne_yksikko) {
    g = ggplot(df1[df1$time_to_dg_mo > -5.1 & df1$tutkimus_lyhenne_yksikko == j,], aes(x = time_to_dg_mo, y = tulos_norm, group=time_to_dg_mo)) +
      geom_jitter(width = 0.2, size = 0.2, color = "black") +
      geom_boxplot(alpha = 0.3, fill = cols1, alpha = 0.3, color = "black", outliers = FALSE) +
      geom_hline(yintercept = median(df1[df1$tutkimus_lyhenne_yksikko == j & df1$time_to_dg_mo == min(df1[df1$time_to_dg_mo > -5.1,]$time_to_dg_mo, na.rm = TRUE),]$tulos_norm, na.rm = TRUE),
                 size = 1, alpha = 1, color = "#e41a1c", lty="11") +
      labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
           y=gsub("^L-", "", j) , x=paste0("Time to diagnosis (years)")) +
      scale_x_continuous(breaks = seq(from = -5, to = -1, by = 1), labels = seq(from = -5, to = -1, by = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title = element_text(size=14, colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size=14, face="bold", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none")
    g
    # Export
    ggsave(plot = g,
           filename = paste0(results, "/", i, "/lab_boxplots/", janitor::make_clean_names(j), "_norm.png"),
           height = 4, width = 4, units = "in", dpi = 300)
  }
}
