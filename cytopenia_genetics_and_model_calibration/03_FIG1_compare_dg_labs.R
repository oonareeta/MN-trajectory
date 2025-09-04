# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Loop over diseases
pvalue_df_labs <- data.frame()
diseases = c("MDS", "MF", "de novo AML")
for (i in diseases) {
  
  # Collect data for disease
  df = fread(paste0(export, "/lagged_data_", i, "3.csv"), nrows = 5) %>%
    dplyr::select(-c(contains("_val_"), contains("ferrit"), contains("_trend_"), contains("l_blast"),
                     contains("l_myelos"), contains("l_meta"), contains("crea"), contains("pt_gf"), contains("p_tt"),
                     contains("rows"), contains("event_1y"), contains("retic")))
  tmp = names(df)
  df = fread(paste0(export, "/lagged_data_", i, "3.csv"), select = tmp)
  df = df %>%
    dplyr::mutate(disease = ifelse(disease == 1, i, "Healthy")) %>%
    mutate(time_to_dg = -time_to_dg) %>%
    dplyr::filter((disease == "Healthy" & time_to_dg<=(-1*365)) | (time_to_dg>=(-2*365) & time_to_dg<(-90)))
  
  
  # Melt
  df = df %>%
    reshape2::melt(id.vars = c("henkilotunnus", "time_to_dg", "disease", "sukupuoli_selite", "age"), variable.name = "tutkimus_lyhenne_yksikko", value.name = "tulos_norm")
  
  # Summarise median by patient and laboratory test
  df = df %>%
    dplyr::group_by(henkilotunnus, tutkimus_lyhenne_yksikko) %>%
    mutate(tulos_norm = median(tulos_norm, na.rm=TRUE)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  
  ############## WILCOXON TESTS ##############
  
  
  ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  multiple_t_tests_p_value <- NULL
  multiple_t_tests_p_value <- lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) wilcox.test(df[df$tutkimus_lyhenne_yksikko == x,]$tulos_norm ~ df[df$tutkimus_lyhenne_yksikko == x,]$disease, na.rm=TRUE))
  ### P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ### Create a matrix and dataframe of the p-values
  pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
    #### Add the p values to a new dataframe
    p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"))
  
  ## Add medians
  pvalue_df1$median0 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) median(df[df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  pvalue_df1$median2 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) median(df[df$disease == i & df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  pvalue_df1$median1 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) median(df[df$disease == "Healthy" & df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  pvalue_df1$mean0 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) mean(df[df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  pvalue_df1$mean2 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) mean(df[df$disease == i & df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  pvalue_df1$mean1 = unlist(lapply(unique(df$tutkimus_lyhenne_yksikko), function(x) mean(df[df$disease == "Healthy" & df$tutkimus_lyhenne_yksikko == x,]$tulos_norm, na.rm=TRUE)))
  
  pvalue_df1 = pvalue_df1 %>%
    ### Rownames to column
    mutate(Lab = unique(df$tutkimus_lyhenne_yksikko),
           Disease = i) %>%
    rename(pvalue = p.value) %>%
    arrange(pvalue)
  
  # Save data
  pvalue_df_labs <- rbind(pvalue_df_labs, pvalue_df1)
  
}

# Save
saveRDS(pvalue_df_labs, paste0(export, "/wilcoxon_labs.rds"))
pvalue_df_labs = readRDS(paste0(export, "/wilcoxon_labs.rds"))


############## PLOTS ##############


# Prepare data for balloonplot
## Calculate FC and -log10 P value
names1 = readxl::read_xlsx(paste0(export, "/lab_test_names.xlsx"))
pvalue_df_labs1 <- pvalue_df_labs %>%
  dplyr::left_join(names1) %>%
  dplyr::select(-Lab) %>%
  dplyr::rename(Lab = Lab1) %>%
  dplyr::filter(!str_detect(Lab, "RDW-SD|L-Blast")) %>%
  mutate(Disease = ifelse(Disease == "MF", "Primary MF", Disease),
         Disease = paste0(toupper(substr(Disease, 1, 1)), substr(Disease, 2, nchar(Disease))),
         Lab = ifelse(Lab=="B-La (mm/h)", "B-ESR (mm/h)", Lab),
         FC = ifelse(median1 == 0 & median2 == 0, 1,
                     ifelse(median2 == 0, 1/median1,
                            ifelse(median1 == 0, median2+1, median2 / median1))),
         p_adjusted = ifelse(p_adjusted == 0, min(p_adjusted[p_adjusted>0], na.rm = TRUE), p_adjusted),
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
                         ifelse(FC<0, log10(-FC), log10(FC))),
         neg_log10_p_adj = ifelse(is.na(FC_log), 0, neg_log10_p_adj))

# Plot
p <- ggballoonplot(pvalue_df_labs1, x = "Lab", y = "Disease",
                   fill = "FC_log",
                   size = "neg_log10_p_adj",
                   # size.range = c(1, 10),
                   ggtheme = theme_bw()) +
  scale_size(
    breaks = c(-log10(0.05), 2, 3),
    range = c(1, 7),
    labels = c("*", "**", "***"),
    limits = c(-log10(0.05), max(pvalue_df_labs1$neg_log10_p_adj, na.rm = TRUE))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0) +
  guides(size = guide_legend(title="Log10 AdjP", nrow = 3, title.vjust = 0.5),
         fill = guide_colorbar(title="Log10 FC", title.vjust = 0.75)) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
        axis.text.x = element_text(colour="black", angle=45, hjust=1),
        axis.title.x = element_text(size=12, colour="black", face="bold"),
        axis.text.y = element_text(colour="black"),
        legend.key.size = unit(0.17, 'in'),
        legend.title = element_text(colour="black", face="bold"),
        legend.position = "right", legend.box="vertical",
        legend.margin=margin()); p
ggsave(plot = p, filename = paste0(results, "/Balloonplot_disease.png"), width = 8, height = (round(length(unique(pvalue_df_labs1$Disease))/5)+2), units = "in", dpi = 300)
