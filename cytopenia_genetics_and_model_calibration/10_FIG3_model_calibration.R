# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Load test data
disease = c("any_MN", "de_novo_AML", "MDS", "MF")

# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv"))
  
  # Select only disease patients
  df_disease = df %>%
    dplyr::rename(
      disease = disease_status) %>%
    dplyr::filter(disease == 1) %>%
    dplyr::mutate(risk_score = log(risk_score, 10),
                  time_to_dg_year = -time_to_dg / 365.24)
  
  # Folder
  if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  # Plot
  g = ggplot(df_disease, aes(x = time_to_dg_year, y = risk_score)) +
    geom_point(size = 0.2, color = "black") +
    geom_smooth(method = "loess") +
    labs(title = paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep=""),
         y="Risk score probability (LOG)",
         x="Time to diagnosis (years)") +
    stat_cor(method = "spearman", label.x = -5,
             label.y = max(df_disease$risk_score, na.rm = TRUE)-0.1) +
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
          legend.position = "none"); g
  
  # Export
  ggsave(plot = g,
         filename = paste0(results, "/", j, "/Calibration_plot_probability_time_to_dg.png"),
         height = 4, width = 4, units = "in", dpi = 300)
}



################ Prediction probability by mutations and karyotype alterations ################


# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = fread(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv"))
  
  # Folder
  if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  if (i == "any_MN") {
    # Load cytogenetics and mutation data
    df_c = fread(paste0(export, "/chromosomes_lab_MDS.csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct() %>%
      full_join(fread(paste0(export, "/chromosomes_lab_de novo AML.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct()) %>%
      full_join(fread(paste0(export, "/chromosomes_lab_MF.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct())
    df_m = fread(paste0(export, "/mutations_lab_MDS.csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct() %>%
      full_join(fread(paste0(export, "/mutations_lab_de novo AML.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct()) %>%
      full_join(fread(paste0(export, "/mutations_lab_MF.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct())
  } else {
    # Load cytogenetics and mutation data
    df_c = fread(paste0(export, "/chromosomes_lab_", j, ".csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct()
    df_m = fread(paste0(export, "/mutations_lab_", j, ".csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct()
  }
  
  # Remove rows where all variables from 2 to ncol are ""
  df_c <- df_c[!apply(df_c[, 2:ncol(df_c)], 1, function(row) all(row == "")), ]
  df_m <- df_m[!apply(df_m[, 2:ncol(df_m)], 1, function(row) all(row == "")), ]
  
  # Names of alterations
  chrm_pos_counts = names(df_c)[2:ncol(df_c)]
  mut_pos_counts = names(df_m)[2:ncol(df_m)]
  
  # Join predictions
  df_c = df_c %>%
    dplyr::inner_join(df %>%
                        dplyr::rename(henkilotunnus = patient_id))
  df_m = df_m %>%
    dplyr::inner_join(df %>%
                        dplyr::rename(henkilotunnus = patient_id))
  
  
  ############## WILCOXON TESTS ##############
  
  
  # Loop karyotypes
  pvalue_df_chrms <- data.frame()
  df_c1 = df_c
  
  ignore_variables1 <- lapply(chrm_pos_counts, function(x) max(df_c1$risk_score[df_c1[[x]] == "NEG"], na.rm=TRUE))
  ignore_variables1 <- chrm_pos_counts[ignore_variables1 == "-Inf"]
  chrm_pos_counts1 <- chrm_pos_counts[!chrm_pos_counts %in% ignore_variables1]
  ignore_variables1 <- lapply(chrm_pos_counts, function(x) max(df_c1$risk_score[df_c1[[x]] == "POS"], na.rm=TRUE))
  ignore_variables1 <- chrm_pos_counts[ignore_variables1 == "-Inf"]
  chrm_pos_counts1 <- chrm_pos_counts1[!chrm_pos_counts1 %in% ignore_variables1]
  
  # Replace "" with NA
  df_c1[df_c1 == ""] <- NA
  
  
  if (length(chrm_pos_counts1) > 0) {
    ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
    multiple_t_tests_p_value <- NULL
    multiple_t_tests_p_value <- lapply(chrm_pos_counts1, function(x) wilcox.test(df_c1$risk_score ~ df_c1[[x]], na.rm=TRUE))
    ### P-values can be extracted from the result object
    pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
    ### Create a matrix and dataframe of the p-values
    pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
      #### Add the p values to a new dataframe
      p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
      #### Add also the t values, 95%CI to the same dataframe
      median0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$risk_score, na.rm=TRUE))),
      median2 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$risk_score[df_c1[[x]] == "POS"], na.rm=TRUE))),
      median1 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$risk_score[df_c1[[x]] == "NEG"], na.rm=TRUE))),
      median2_0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$risk_score[df_c1[[x]] == "POS"], na.rm=TRUE))),
      median1_0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$risk_score[df_c1[[x]] == "NEG"], na.rm=TRUE))),
      n_pos = lapply(chrm_pos_counts1, function(x) sum(df_c1[[x]] == "POS", na.rm=TRUE) ),
      n_neg = lapply(chrm_pos_counts1, function(x) sum(df_c1[[x]] == "NEG", na.rm=TRUE) )
    ) %>%
      mutate(disease = i,
             gene = chrm_pos_counts1) %>%
      rename(pvalue = p.value) %>%
      arrange(pvalue)
    rownames(pvalue_df1) = NULL
    
    # Save data
    pvalue_df_chrms <- rbind(pvalue_df_chrms, pvalue_df1)
    
  }
  
  
  # Loop mutations
  pvalue_df_muts <- data.frame()
  df_m1 = df_m
  
  # Ignore variables with all NEGs
  ignore_variables1 <- lapply(mut_pos_counts, function(x) max(df_m1$risk_score[df_m1[[x]] == "NEG"], na.rm=TRUE))
  ignore_variables1 <- mut_pos_counts[ignore_variables1 == "-Inf"]
  mut_pos_counts1 <- mut_pos_counts[!mut_pos_counts %in% ignore_variables1]
  ignore_variables1 <- lapply(mut_pos_counts, function(x) max(df_m1$risk_score[df_m1[[x]] == "POS"], na.rm=TRUE))
  ignore_variables1 <- mut_pos_counts[ignore_variables1 == "-Inf"]
  mut_pos_counts1 <- mut_pos_counts1[!mut_pos_counts1 %in% ignore_variables1]
  
  # Replace "" with NA
  df_m1[df_m1 == ""] <- NA
  
  if (length(mut_pos_counts1) > 0) {
    ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
    multiple_t_tests_p_value <- NULL
    multiple_t_tests_p_value <- lapply(mut_pos_counts1, function(x) wilcox.test(df_m1$risk_score ~ df_m1[[x]], na.rm=TRUE))
    ### P-values can be extracted from the result object
    pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
    ### Create a matrix and dataframe of the p-values
    pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
      #### Add the p values to a new dataframe
      p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
      #### Add also the t values, 95%CI to the same dataframe
      median0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$risk_score, na.rm=TRUE))),
      median2 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$risk_score[df_m1[[x]] == "POS"], na.rm=TRUE))),
      median1 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$risk_score[df_m1[[x]] == "NEG"], na.rm=TRUE))),
      median2_0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$risk_score[df_m1[[x]] == "POS"], na.rm=TRUE))),
      median1_0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$risk_score[df_m1[[x]] == "NEG"], na.rm=TRUE))),
      n_pos = lapply(mut_pos_counts1, function(x) sum(df_m1[[x]] == "POS", na.rm=TRUE) ),
      n_neg = lapply(mut_pos_counts1, function(x) sum(df_m1[[x]] == "NEG", na.rm=TRUE) )
    ) %>%
      mutate(disease = i,
             gene = mut_pos_counts1) %>%
      rename(pvalue = p.value) %>%
      arrange(pvalue)
    rownames(pvalue_df1) = NULL
    
    # Save data
    pvalue_df_muts <- rbind(pvalue_df_muts, pvalue_df1)
    
  }
  
  # Save
  saveRDS(pvalue_df_chrms, paste0(export, "/wilcoxon_chromosomes_", i, "_model_predictions.rds"))
  saveRDS(pvalue_df_muts, paste0(export, "/wilcoxon_mutations_", i, "_model_predictions.rds"))
  
}


############## PLOTS ##############


for (i in disease) {
  
  
  # Folder
  if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  # Read data
  pvalue_df_chrms = readRDS(paste0(export, "/wilcoxon_chromosomes_", i, "_model_predictions.rds"))
  if (nrow(pvalue_df_chrms) > 0) {
    pvalue_df_chrms = pvalue_df_chrms %>%
      # Rename alterations
      dplyr::mutate(
        gene = gsub("trisomy", "+",
                    gsub("_", ",",
                         gsub("-y", "-Y",
                              gsub("(monosomy|monosomy_)", "-",
                                   gsub("_kar", " kar",
                                        gsub("karyotype_HR", "High-risk karyotype",
                                             gsub("translocation", "t(",
                                                  gsub("t_", "t(",
                                                       gsub("inv_", "inv(",
                                                            gsub("none", "Normal karyotype",
                                                                 gsub("del_", "del(", gene))))))))))),
        gene = ifelse(str_detect(gene, "(^(t\\(|inv|del))"), paste0(gene, ")"), gene),
        gene = gsub("mono", "Mono",
                    gsub("com", "Com",
                         gsub("t\\(3\\)", "t(3,)",
                              gsub("t\\(11\\)", "t(11,)", gene))))
      )
  }
  
  pvalue_df_muts = readRDS(paste0(export, "/wilcoxon_mutations_", i, "_model_predictions.rds"))
  if (nrow(pvalue_df_muts) > 0) {
    pvalue_df_muts = pvalue_df_muts %>%
      dplyr::mutate(gene = ifelse(gene == "None", "No mutation", gsub("FLT3_", "FLT3-", gene)))
  }
  
  pvalue_df_chrms = rbind(pvalue_df_chrms, pvalue_df_muts)
  
  
  # CHROMOSOMES
  # Prepare data for balloonplot
  ## Calculate FC and -log10 P value
  pvalue_df_chrm3 <- pvalue_df_chrms %>%
    mutate(FC = ifelse(median1_0 == 0 & median2_0 == 0, 1,
                       ifelse(median2_0 == 0, 1/median1_0+1,
                              ifelse(median1_0 == 0, median2_0+1, median2_0/median1_0)))) %>%
    mutate(
      neg_log10_p = -log10(pvalue),
      neg_log10_p_adj = -log10(p_adjusted),
      p_adjusted_cat = ifelse(p_adjusted<0.001, 0.001,
                              ifelse(p_adjusted<0.01, 0.01,
                                     ifelse(p_adjusted<0.05, 0.05,
                                            ifelse(p_adjusted<0.1, 0.1, "ns")))),
      pvalue_cat = ifelse(pvalue<0.001, 0.001,
                          ifelse(pvalue<0.01, 0.01,
                                 ifelse(pvalue<0.05, 0.05,
                                        ifelse(pvalue<0.1, 0.1, "ns")))),
      FC_log = ifelse(FC == 0, NA,
                      ifelse(FC<0, log(-FC), log(FC))),
      neg_log10_p_adj = ifelse(is.na(FC_log), 0, neg_log10_p_adj))
  
  # Order by value
  pvalue_df_chrm2 = pvalue_df_chrm3 %>%
    dplyr::mutate(FC1 = ifelse(FC < 1, p_adjusted, 1/p_adjusted)) %>%
    dplyr::arrange(FC1) %>%
    dplyr::select(gene); ordered_columns2 = pvalue_df_chrm2$gene
  
  # Replace values lower than -3 with -3 and higher than 3 with 3
  pvalue_df_chrm3 = pvalue_df_chrm3 %>%
    dplyr::mutate(
      FC_log = ifelse(FC_log < -1.5, -1.5, ifelse(FC_log > 1.5, 1.5, FC_log)),
      disease = "")
  
  # Color genes
  cols1 = ifelse(ordered_columns2 %in% unique(pvalue_df_muts$gene), "#a6761d", "black")
  
  
  # Balloonplot
  p <- ggballoonplot(pvalue_df_chrm3, x = "disease", y = "gene",
                     fill = "FC_log",
                     size = "neg_log10_p_adj",
                     ggtheme = theme_bw()) +
    scale_size(breaks = c(-log10(0.05), -log10(0.01), -log10(0.001)), #c(1, 2, 3),
               labels = c("<0.05", "<0.01", "<0.001"),
               range = c(2, 10),
               limits = c(-log10(0.05), max(pvalue_df_chrm3$neg_log10_p_adj, na.rm = TRUE))) +
    scale_y_discrete(limits = ordered_columns2) +
    xlab(paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep="")) +
    ylab("") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         breaks = c(-2, -1, 0, 1, 2),
                         midpoint = 0) +
    guides(size = guide_legend(title="LOG10 AdjP", nrow = 1, title.vjust = 0.5),
           fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
    theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
          axis.text.x = element_text(colour="black", angle=45, hjust=1),
          axis.title.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(colour=cols1),
          panel.background = element_rect(fill='transparent', color = NA),
          rect = element_rect(fill = "transparent"),
          legend.position = "bottom", legend.box="vertical", legend.margin=margin(),
          plot.margin = unit(c(1,3.5,1,2.5), "cm")); p
  
  # Balloonplot
  if (nrow(pvalue_df_chrm3)==1) {
    p <- ggballoonplot(pvalue_df_chrm3, x = "disease", y = "gene",
                       fill = "FC_log",
                       size = "neg_log10_p_adj",
                       ggtheme = theme_bw()) +
      scale_size(breaks = c(-log10(0.05), -log10(0.01), -log10(0.001)), #c(1, 2, 3),
                 labels = c("<0.05", "<0.01", "<0.001"),
                 range = c(2, 10),
                 limits = c(-log10(0.05), max(pvalue_df_chrm3$neg_log10_p_adj, na.rm = TRUE))) +
      scale_y_discrete(limits = ordered_columns2) +
      xlab(paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep="")) +
      ylab("") +
      guides(size = guide_legend(title="LOG10 AdjP", nrow = 1, title.vjust = 0.5),
             fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
      theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
            axis.text.x = element_text(colour="black", angle=45, hjust=1),
            axis.title.x = element_text(size=12, colour="black", face="bold"),
            axis.text.y = element_text(colour=cols1),
            panel.background = element_rect(fill='transparent', color = NA),
            rect = element_rect(fill = "transparent"),
            legend.position = "bottom", legend.box="vertical", legend.margin=margin(),
            plot.margin = unit(c(1,3.5,1,2.5), "cm")); p
  }
  
  ggsave(plot = p, filename = paste0(results, "/", j, "/Balloonplot_chromosomes_model_prediction.png"), bg = "white",
         width = 4.5, height = (floor(length(unique(pvalue_df_chrm3$gene))/5)+3.5), units = "in", dpi = 300)
  
}
