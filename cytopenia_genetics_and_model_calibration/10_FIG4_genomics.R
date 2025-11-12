# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")

myhclust <- function(x){
  cor <- cor(as.matrix(x), method="spearman", use="pairwise.complete.obs")
  cor2 <- as.dist(1-cor)
  hclust(d=cor2, method="ward.D2")
}

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# Loop over diseases
diseases = c("MDS", "MF", "de novo AML") #c("AML", "MDS", "MF", "de novo AML", "secondary AML")
for (i in diseases) {
  
  # Collect data for disease
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read lab data
  df = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::filter(time_to_dg<-90 & time_to_dg>=-730)
  
  # Summarise median by patient and laboratory test
  df = df %>%
    dplyr::group_by(henkilotunnus, tutkimus_lyhenne_yksikko) %>%
    mutate(tulos_norm = mean(tulos_norm, na.rm=TRUE)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  
  # Molecular genetics
  if (i %in% c("AML", "de novo AML", "secondary AML")) {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/AML/ELN.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), translocation3, translocation11, complex_karyotype, monosomal_karyotype, karyotype_HR) %>%
      dplyr::mutate_at(.vars = c("translocation3", "translocation11", "complex_karyotype", "monosomal_karyotype", "karyotype_HR"), function(x) ifelse(is.na(x), NA,
                                                                                                                                                      ifelse(x==1, "POS", "NEG")))
    
  } else if (i %in% c("MDS")) {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/MDS/IPSSR.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), translocation3, translocation11, complex_karyotype, monosomal_karyotype, karyotype_HR) %>%
      dplyr::mutate_at(.vars = c("translocation3", "translocation11", "complex_karyotype", "monosomal_karyotype", "karyotype_HR"), function(x) ifelse(is.na(x), NA,
                                                                                                                                                      ifelse(x==1, "POS", "NEG")))
    
  } else {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/MF/ELN.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), complex_karyotype, monosomal_karyotype) %>%
      dplyr::mutate_at(.vars = c("complex_karyotype", "monosomal_karyotype"), function(x) ifelse(is.na(x), NA,
                                                                                                 ifelse(x==1, "POS", "NEG")))
  }
  
  # Edit none
  mut = mut %>%
    dplyr::mutate(None = ifelse(is.na(None) & !is.na(ASXL1), "NEG", None))
  
  # Mutation data
  muts = mut %>%
    dplyr::select(henkilotunnus, which(names(mut) == "ABL1") : which(names(mut) == "ZRSR2"))
  
  
  ## Keep columns with ≥10 POS
  # Function to count POS per column
  count_pos_per_column <- function(df) {
    sapply(df, function(x) sum(x == "POS", na.rm=TRUE))
  }
  
  # Count NAs in the sample data frame
  n1 = pmin(15, 0.05*nrow(muts))
  
  
  ## Chromosomes
  chrm_pos_counts = count_pos_per_column(chrms)
  chrm_pos_counts = chrm_pos_counts[chrm_pos_counts>=n1]
  chrm_pos_counts = chrm_pos_counts[names(chrm_pos_counts)[grep(pattern = "der_", x = names(chrm_pos_counts), invert = TRUE)]]
  ## Mutations
  mut_pos_counts = count_pos_per_column(muts)
  mut_pos_counts = mut_pos_counts[mut_pos_counts>=n1]
  
  # Combine
  df_c = df %>%
    dplyr::inner_join(chrms %>% dplyr::select(henkilotunnus, all_of(names(chrm_pos_counts))))
  
  df_c = df_c %>%
    dplyr::select(-contains("t_11q_3"), -starts_with("der"))
  df_m = df %>%
    dplyr::inner_join(muts %>% dplyr::select(henkilotunnus, all_of(names(mut_pos_counts))))
  
  # Save
  fwrite(df_m, paste0(export, "/mutations_lab_", i, ".csv"))
  fwrite(df_c, paste0(export, "/chromosomes_lab_", i, ".csv"))
  
  
  ############## WILCOXON TESTS ##############
  
  
  # Loop karyotypes
  pvalue_df_chrms <- data.frame()
  for (lab1 in 1:length(unique(df$tutkimus_lyhenne_yksikko))) {
    print(paste0(lab1, "/", length(unique(df$tutkimus_lyhenne_yksikko))))
    df_c1 = df_c %>%
      dplyr::filter(tutkimus_lyhenne_yksikko == unique(df$tutkimus_lyhenne_yksikko)[lab1]) %>%
      dplyr::mutate(tulos_norm_scaled = scales::rescale(tulos_norm, to = c(0,1)))
    
    # Ignore variables with all NEGs
    ignore_variables1 <- lapply(names(chrm_pos_counts), function(x) max(df_c1$tulos_norm[df_c1[[x]] == "NEG"], na.rm=TRUE))
    ignore_variables1 <- names(chrm_pos_counts)[ignore_variables1 == "-Inf"]
    chrm_pos_counts1 <- names(chrm_pos_counts)[!names(chrm_pos_counts) %in% ignore_variables1]
    ignore_variables1 <- lapply(names(chrm_pos_counts), function(x) max(df_c1$tulos_norm[df_c1[[x]] == "POS"], na.rm=TRUE))
    ignore_variables1 <- names(chrm_pos_counts)[ignore_variables1 == "-Inf"]
    chrm_pos_counts1 <- chrm_pos_counts1[!chrm_pos_counts1 %in% ignore_variables1]
    
    
    if (length(chrm_pos_counts1) > 0) {
      ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
      multiple_t_tests_p_value <- NULL
      multiple_t_tests_p_value <- lapply(chrm_pos_counts1, function(x) wilcox.test(df_c1$tulos_norm ~ df_c1[[x]], na.rm=TRUE))
      ### P-values can be extracted from the result object
      pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
      ### Create a matrix and dataframe of the p-values
      pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
        #### Add the p values to a new dataframe
        p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
        #### Add also the t values, 95%CI to the same dataframe
        median0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$tulos_norm, na.rm=TRUE))),
        median2 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$tulos_norm[df_c1[[x]] == "POS"], na.rm=TRUE))),
        median1 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$tulos_norm[df_c1[[x]] == "NEG"], na.rm=TRUE))),
        median2_0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$tulos_norm_scaled[df_c1[[x]] == "POS"], na.rm=TRUE))),
        median1_0 = unlist(lapply(chrm_pos_counts1, function(x) median(df_c1$tulos_norm_scaled[df_c1[[x]] == "NEG"], na.rm=TRUE))),
        n_pos = lapply(chrm_pos_counts1, function(x) sum(df_c1[[x]] == "POS", na.rm=TRUE) ),
        n_neg = lapply(chrm_pos_counts1, function(x) sum(df_c1[[x]] == "NEG", na.rm=TRUE) )
      ) %>%
        ### Rownames to column
        rownames_to_column() %>%
        mutate(rowname = unique(df$tutkimus_lyhenne_yksikko)[lab1],
               gene = chrm_pos_counts1) %>%
        rename(Lab = rowname,
               pvalue = p.value) %>%
        arrange(pvalue)
      rownames(pvalue_df1) = NULL
      
      # Add ignored variables as NA
      for (h in ignore_variables1) {
        pvalue_df2 = pvalue_df1[1,]
        pvalue_df2[1, ] <- NA
        pvalue_df2 = pvalue_df2 %>%
          dplyr::mutate(gene = h,
                        Lab = unique(df$tutkimus_lyhenne_yksikko)[lab1])
        pvalue_df1 <- rbind(pvalue_df1, pvalue_df2)
      }
      
      # Save data
      pvalue_df_chrms <- rbind(pvalue_df_chrms, pvalue_df1)
    }
    
  }
  
  # Loop mutations
  pvalue_df_muts <- data.frame()
  for (lab1 in 1:length(unique(df$tutkimus_lyhenne_yksikko))) {
    print(paste0(lab1, "/", length(unique(df$tutkimus_lyhenne_yksikko))))
    df_m1 = df_m %>%
      dplyr::filter(tutkimus_lyhenne_yksikko == unique(df$tutkimus_lyhenne_yksikko)[lab1]) %>%
      dplyr::mutate(tulos_norm_scaled = scales::rescale(tulos_norm, to = c(0,1)))
    
    # Ignore variables with all NEGs
    ignore_variables1 <- lapply(names(mut_pos_counts), function(x) max(df_m1$tulos_norm[df_m1[[x]] == "NEG"], na.rm=TRUE))
    ignore_variables1 <- names(mut_pos_counts)[ignore_variables1 == "-Inf"]
    mut_pos_counts1 <- names(mut_pos_counts)[!names(mut_pos_counts) %in% ignore_variables1]
    ignore_variables1 <- lapply(names(mut_pos_counts), function(x) max(df_m1$tulos_norm[df_m1[[x]] == "POS"], na.rm=TRUE))
    ignore_variables1 <- names(mut_pos_counts)[ignore_variables1 == "-Inf"]
    mut_pos_counts1 <- mut_pos_counts1[!mut_pos_counts1 %in% ignore_variables1]
    
    if (length(mut_pos_counts1) > 0) {
      ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
      multiple_t_tests_p_value <- NULL
      multiple_t_tests_p_value <- lapply(mut_pos_counts1, function(x) wilcox.test(df_m1$tulos_norm ~ df_m1[[x]], na.rm=TRUE))
      ### P-values can be extracted from the result object
      pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
      ### Create a matrix and dataframe of the p-values
      pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
        #### Add the p values to a new dataframe
        p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
        #### Add also the t values, 95%CI to the same dataframe
        median0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$tulos_norm, na.rm=TRUE))),
        median2 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$tulos_norm[df_m1[[x]] == "POS"], na.rm=TRUE))),
        median1 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$tulos_norm[df_m1[[x]] == "NEG"], na.rm=TRUE))),
        median2_0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$tulos_norm_scaled[df_m1[[x]] == "POS"], na.rm=TRUE))),
        median1_0 = unlist(lapply(mut_pos_counts1, function(x) median(df_m1$tulos_norm_scaled[df_m1[[x]] == "NEG"], na.rm=TRUE))),
        n_pos = lapply(mut_pos_counts1, function(x) sum(df_m1[[x]] == "POS", na.rm=TRUE) ),
        n_neg = lapply(mut_pos_counts1, function(x) sum(df_m1[[x]] == "NEG", na.rm=TRUE) )
      ) %>%
        ### Rownames to column
        rownames_to_column() %>%
        mutate(rowname = unique(df$tutkimus_lyhenne_yksikko)[lab1],
               gene = mut_pos_counts1) %>%
        rename(Lab = rowname,
               pvalue = p.value) %>%
        arrange(pvalue)
      rownames(pvalue_df1) = NULL
      
      # Add ignored variables as NA
      for (h in ignore_variables1) {
        pvalue_df2 = pvalue_df1[1,]
        pvalue_df2[1, ] <- NA
        pvalue_df2 = pvalue_df2 %>%
          dplyr::mutate(gene = h,
                        Lab = unique(df$tutkimus_lyhenne_yksikko)[lab1])
        pvalue_df1 <- rbind(pvalue_df1, pvalue_df2)
      }
      
      # Save data
      pvalue_df_muts <- rbind(pvalue_df_muts, pvalue_df1)
    }
    
  }
  
  # Save
  saveRDS(pvalue_df_chrms, paste0(export, "/wilcoxon_chromosomes_", i, ".rds"))
  saveRDS(pvalue_df_muts, paste0(export, "/wilcoxon_mutations_", i, ".rds"))
  
  
  ############## PLOTS ##############
  
  
  # Read data
  pvalue_df_chrms = readRDS(paste0(export, "/wilcoxon_chromosomes_", i, ".rds")) %>%
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
  
  pvalue_df_muts = readRDS(paste0(export, "/wilcoxon_mutations_", i, ".rds")) %>%
    dplyr::mutate(gene = ifelse(gene == "None", "No mutation", gsub("FLT3_", "FLT3-", gene)))
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
      neg_log10_p_adj = -log10(pvalue),
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
  
  
  # Hclusters
  ## CHROMOSOMES
  # Calculate the distance matrix using Euclidean distance
  pvalue_df_chrm1 = dcast(pvalue_df_chrm3 %>%
                            dplyr::mutate(
                              FC_log = ifelse(FC_log < -2, -2, ifelse(FC_log > 2, 2, FC_log)),
                              neg_log10_p_adj = ifelse(neg_log10_p_adj > 5, 5, neg_log10_p_adj),
                              FC = neg_log10_p_adj*FC_log) %>%
                            dplyr::select(Lab, gene, FC), gene~Lab, value.var = "FC")
  dist_matrix <- dist(pvalue_df_chrm1[,2:ncol(pvalue_df_chrm1)], method = "euclidean")
  # Perform hierarchical clustering using the ward.D2 method
  hclust_result <- myhclust(dist_matrix)
  # Get the order of columns based on the hierarchical clustering
  ordered_columns1 <- hclust_result$order
  cols1 = ifelse(pvalue_df_chrm1$gene[ordered_columns1] %in% unique(pvalue_df_muts$gene), "black", "#a6761d")
  
  if (length(ordered_columns1) > 2) {
    # Plot clusters
    png(filename = paste0(results, "/", i, "/Balloonplot_chromosomes_clusters.png"), width = 6, height = 6, units = "in", res = 300)
    plot(hclust_result)
    dev.off()
  }
  
  ## LABS
  # Compute the Euclidean distance between columns
  pvalue_df_chrm2 = dcast(pvalue_df_chrm3 %>%
                            dplyr::mutate(
                              FC_log = ifelse(FC_log < -2, -2, ifelse(FC_log > 2, 2, FC_log)),
                              neg_log10_p_adj = ifelse(neg_log10_p_adj > 5, 5, neg_log10_p_adj),
                              FC = neg_log10_p_adj*FC_log) %>%
                            dplyr::select(Lab, gene, FC), Lab~gene, value.var = "FC")
  dist_matrix <- dist(pvalue_df_chrm2[,2:ncol(pvalue_df_chrm2)], method = "euclidean")
  # Perform hierarchical clustering using the ward.D2 method
  hclust_result <- myhclust(dist_matrix)
  # Get the order of columns based on the hierarchical clustering
  ordered_columns2 <- hclust_result$order
  
  # Replace values lower than -3 with -3 and higher than 3 with 3
  pvalue_df_chrm3 = pvalue_df_chrm3 %>%
    dplyr::mutate(
      FC_log = ifelse(FC_log < -1.5, -1.5, ifelse(FC_log > 1.5, 1.5, FC_log)))
  
  
  # Balloonplot
  p <- ggballoonplot(pvalue_df_chrm3, x = "Lab", y = "gene",
                     fill = "FC_log",
                     size = "neg_log10_p_adj",
                     ggtheme = theme_bw()) +
    scale_size(breaks = c(-log10(0.05), -log10(0.01), -log10(0.001)), #c(1, 2, 3),
               labels = c("<0.05", "<0.01", "<0.001"),
               range = c(2, 10),
               limits = c(-log10(0.05), max(pvalue_df_chrm3$neg_log10_p_adj, na.rm = TRUE))) +
    scale_x_discrete(limits = names(pvalue_df_chrm1)[2:ncol(pvalue_df_chrm1)][ordered_columns2]) +
    scale_y_discrete(limits = names(pvalue_df_chrm2)[2:ncol(pvalue_df_chrm2)][ordered_columns1]) +
    ylab("") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         breaks = c(-2, -1, 0, 1, 2),
                         midpoint = 0) +
    guides(size = guide_legend(title="LOG10 P", nrow = 1, title.vjust = 0.5),
           fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
    theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
          axis.text.x = element_text(colour="black", angle=45, hjust=1),
          axis.title.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(colour=cols1),
          panel.background = element_rect(fill='transparent', color = NA),
          rect = element_rect(fill = "transparent"),
          legend.position = "bottom", legend.box="horizontal", legend.margin=margin()); p
  ggsave(plot = p, bg = "white", filename = paste0(results, "/", i, "/Balloonplot_chromosomes.png"), width = 8, height = (floor(length(unique(pvalue_df_chrm3$gene))/5)+3.5), units = "in", dpi = 300)

  
}

