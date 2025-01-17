# Link cytogenetics and genomics data to temporal laboratory values

# Author: Oscar Brück

# Libraries
source("./library.R")

# General parameters
source("./parameters")


############################################################


# Load test data
disease = c("any_MN", "de_novo_AML", "MDS", "primary_MF")

"mounts/research/husdatalake/disease/processed_data/Preleukemia/de novo AML_lab_demo.rds"
# Prediction probability by mutations and karyotype alterations
# Loop over diseases
for (i in disease) {
  i = "any_MN"
  # Load test data
  df_train = fread(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona2/results/final_model/", i, "_train_data_with_final_model_predictions.csv"), select = c("henkilotunnus", "time_to_dg", "disease", "b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm", "b_wbc_e9_l_tulos_norm"))
  df_val = fread(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona2/results/final_model/", i, "_validation_data_with_final_model_predictions.csv"), select = c("henkilotunnus", "time_to_dg", "disease", "b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm", "b_wbc_e9_l_tulos_norm"))
  df_test = fread(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona2/results/final_model/", i, "_test_data_with_final_model_predictions.csv"), select = c("henkilotunnus", "time_to_dg", "disease", "b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm", "b_wbc_e9_l_tulos_norm"))
  df = rbind(df_train, df_val, df_test)
  
  # Folder
  if (i == "primary_MF") {
    j = "MF"
  } else if (i == "de_novo_AML") {
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
  alterations = c(chrm_pos_counts, mut_pos_counts)
  labs = c("b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm", "b_wbc_e9_l_tulos_norm")
  
  # # Add healthy
  # df_m_template = df_m[1,2:ncol(df_m)]
  # df_m_template[1,2:ncol(df_m_template)] = "NEG"
  # healthy_m = df %>%
  #   dplyr::filter(disease == 0) %>%
  #   distinct(henkilotunnus) %>%
  #   bind_cols(df_m_template)
  # df_c_template = df_c[1,2:ncol(df_c)]
  # df_c_template[1,2:ncol(df_c_template)] = "NEG"
  # healthy_c = df %>%
  #   dplyr::filter(disease == 0) %>%
  #   distinct(henkilotunnus) %>%
  #   bind_cols(df_c_template)
  # 
  # df_m = df_m %>%
  #   rbind(healthy_m)
  # df_c = df_c %>%
  #   rbind(healthy_c)
  
  # Join predictions
  df_c = df_c %>%
    dplyr::inner_join(df)
  df_m = df_m %>%
    dplyr::inner_join(df)
  
  df1 = full_join(df_m, df_c)
  
  
  df2 = df_m %>%
    # dplyr::filter(disease == 1) %>%
    mutate(time_to_dg = -round(time_to_dg / 365.24)) %>%
    dplyr::select(-disease) %>%
    reshape2::melt(id.vars = c("henkilotunnus", labs, "time_to_dg"), variable.name = "Alteration") %>%
    dplyr::filter(!value == "") %>%
    reshape2::melt(id.vars = c("henkilotunnus", "time_to_dg", "Alteration", "value"), variable.name = "Labs", value.name = "Lab_value") %>%
    distinct()
  df2_1 = df2 %>%
    dplyr::filter(value == "POS") %>%
    group_by(henkilotunnus, time_to_dg, Alteration, Labs) %>%
    summarise(Lab_value = median(Lab_value, na.rm=TRUE)) %>%
    ungroup() %>%
    dplyr::select(-c(henkilotunnus))
  df2_2 = df2 %>%
    dplyr::filter(value == "NEG") %>%
    group_by(henkilotunnus, time_to_dg, Alteration, Labs) %>%
    summarise(Lab_value = median(Lab_value, na.rm=TRUE)) %>%
    ungroup() %>%
    dplyr::select(-c(henkilotunnus)) %>%
    dplyr::rename(Lab_value_neg = Lab_value)
  df2 = full_join(df2_1, df2_2)
  df2 = df2 %>%
    dplyr::filter(!(is.na(Lab_value) | is.na(Lab_value_neg))) %>%
    dplyr::mutate(Lab_value = Lab_value - Lab_value_neg) %>%
    dplyr::select(-Lab_value_neg)
  
  library(rstatix)
  
  # cor_results = df2 %>%
  #   arrange(desc(time_to_dg)) %>%
  #   group_by(Labs) %>%
  #   dplyr::slice(1) %>%
  #   ungroup()
  # cor_results = cor_results %>%
  #   dplyr::filter(Alteration %in% names(table(cor_results$Alteration, cor_results$Labs)[,1][table(cor_results$Alteration, cor_results$Labs)[,1] > 4]))
  cor_results = df2 %>%
    ungroup() %>%
    # dplyr::filter(Alteration %in% cor_results$Alteration) %>%
    group_by(Alteration, Labs) %>%
    rstatix::cor_test(Lab_value, time_to_dg, method = "spearman") %>%
    rstatix::adjust_pvalue(method = "BH")
  cor_results1 = cor_results %>%
    dplyr::filter(p.adj < 0.05 & abs(cor) > 0.05)
  
  df3 = df2 %>%
    dplyr::inner_join(cor_results1)
  
  df3$y_lab = gsub("b ", "B-", gsub("e9 l", "(E9/l)", gsub("_", " ", gsub("_tulos_norm", "", df3$Labs))))
  df3$y_lab = paste0(toupper(substr(df3$y_lab, 1, 3)), substr(df3$y_lab, 4, nchar(df3$y_lab)), ", normalized", sep="")
  
  # CHROMOSOMES
  # Prepare data for balloonplot
  ## Calculate FC and -log10 P value
  df3 <- df3 %>%
    mutate(
      neg_log10_p = -log10(p),
      neg_log10_p_adj = -log10(p.adj),
      neg_log10_p_adj = ifelse(neg_log10_p_adj > 30, 10, neg_log10_p_adj))
  
  myhclust <- function(x){
    cor <- cor(as.matrix(x), method="spearman", use="pairwise.complete.obs")
    cor2 <- as.dist(1-cor)
    hclust(d=cor2, method="ward.D2")
  }
  
  
  
  df4 = df3 %>%
    dplyr::distinct(Labs, Alteration, cor, neg_log10_p_adj)
  
  # Hclusters
  ## CHROMOSOMES
  # Calculate the distance matrix using Euclidean distance
  pvalue_df_chrm1 = dcast(df4 %>%
                            dplyr::select(Labs, Alteration, cor), Alteration~Labs, value.var = "cor")
  dist_matrix <- dist(pvalue_df_chrm1[,2:ncol(pvalue_df_chrm1)], method = "euclidean")
  # Perform hierarchical clustering using the ward.D2 method
  hclust_result <- myhclust(dist_matrix)
  # Get the order of columns based on the hierarchical clustering
  ordered_columns1 <- hclust_result$order
  cols1 = ifelse(pvalue_df_chrm1$Alteration[ordered_columns1] %in% unique(df4$Alteration), "black", "#a6761d")
  
  ## LABS
  # Compute the Euclidean distance between columns
  pvalue_df_chrm2 = dcast(df4 %>%
                            dplyr::select(Labs, Alteration, cor), Labs~Alteration, value.var = "cor")
  # pvalue_df_chrm2[sapply(pvalue_df_chrm2, is.infinite)] <- 0
  dist_matrix <- dist(pvalue_df_chrm2[,2:ncol(pvalue_df_chrm2)], method = "euclidean")
  # Perform hierarchical clustering using the ward.D2 method
  hclust_result <- myhclust(dist_matrix)
  # Get the order of columns based on the hierarchical clustering
  ordered_columns2 <- hclust_result$order
  
  # Balloonplot
  p <- ggballoonplot(df4, x = "Labs", y = "Alteration",
                     fill = "cor",
                     size = "neg_log10_p_adj",
                     # size.range = c(1, 10),
                     ggtheme = theme_bw()) +
    scale_size(breaks = c(-log10(0.05), -log10(0.01), -log10(0.001)), #c(1, 2, 3),
               labels = c("<0.05", "<0.01", "<0.001"),
               range = c(2, 10),
               limits = c(-log10(0.05), max(df4$neg_log10_p_adj, na.rm = TRUE))) +
    scale_y_discrete(limits = pvalue_df_chrm1[,1][ordered_columns1]) +
    scale_x_discrete(limits = pvalue_df_chrm2[,1][ordered_columns2]) +
    xlab(paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep="")) +
    ylab("") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         # breaks = c(-2, -1, 0, 1, 2),
                         midpoint = 0) +
    guides(size = guide_legend(title="LOG10 AdjP", nrow = 1, title.vjust = 0.5),
           fill = guide_colorbar(title="Cor", title.vjust = 0.75)) +
    # font("xy.text", size = 10, color = "black", face="plain") +
    theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
          axis.text.x = element_text(colour="black", angle=45, hjust=1),
          axis.title.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(colour="black"),
          # axis.text.y = element_text(colour=cols1),
          panel.background = element_rect(fill='transparent', color = NA),
          rect = element_rect(fill = "transparent"),
          legend.position = "bottom", legend.box="vertical", legend.margin=margin(),
          plot.margin = unit(c(1,3.5,1,2.5), "cm")); p
  ggsave(plot = p, filename = paste0(results, "/", j, "/Balloonplot_chromosomes_model_prediction1.png"), bg = "white",
         width = 4.5, height = (floor(length(unique(df3$gene))/5)+3.5), units = "in", dpi = 300)
  
  # # Plot
  # g = ggplot(df3, aes(x = time_to_dg, y = log(Lab_value), color = Alteration)) +
  #   geom_smooth(method = "lm", se = FALSE) + #method = "loess",
  #   labs(title = paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep=""),
  #        y="Blood counts",
  #        x="Time to diagnosis (years)") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(size=12, colour = "black"),
  #         axis.text.y = element_text(size=12, colour = "black"),
  #         axis.title = element_text(size=14, colour = "black"),
  #         axis.line = element_line(colour = "black"),
  #         plot.title = element_text(size=14, face="bold", colour = "black"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_blank(),
  #         panel.background = element_blank(),
  #         legend.position = "bottom") +
  #   facet_wrap(facets = df3$Labs); g
  # 
  # # Export
  # ggsave(plot = g,
  #        filename = paste0(results, "/", j, "/Calibration_plot_probability_time_to_dg1.png"),
  #        height = 4, width = 4, units = "in", dpi = 300)
}
