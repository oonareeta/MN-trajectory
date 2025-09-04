# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Loop over diseases
pvalue_df_labs <- data.frame()
diseases = c("any_MN", "MDS", "MF", "de_novo_AML")
df = data.frame()
for (i in diseases) {
  
  # Collect data for disease
  df1 = fread(paste0(scripts, "/oona_new/explainability/comparison_dataframes/youden/", i, "_TP_FN_comparison_labs.csv")) %>%
    as.data.frame() %>%
    dplyr::rename("Variable" = "V1",
                  "N_TP" = "N TP",
                  "N_FN" = "N FN",
                  "TP_median" = "TP median",
                  "FN_median" = "FN median") %>%
    dplyr::mutate(Variable = ifelse(Variable=="B-La (mm/h)", "B-ESR (mm/h)", Variable),
                  p_adjusted = p.adjust(pvalue, "BH"),
                  Disease = gsub("_", " ", i))

  # Rbind
  df = rbind(df, df1)

}

# Export
fwrite(df, paste0(export, "/TP_FN_comparison_labs_combined.csv"))


############## PLOTS ##############


# Prepare data for balloonplot
## Calculate FC and -log10 P value
pvalue_df_labs1 <- df %>%
  dplyr::filter(!str_detect(Variable, "E-Retic|RDW-SD|rows_in_last_month")) %>%
  mutate(Variable = ifelse(Variable == "Age", "age", Variable),
         Variable = ifelse(Variable == "Sex", "sukupuoli_selite", Variable),
         FC = ifelse(TP_median == 0 & FN_median == 0, 1,
                     ifelse(FN_median == 0, TP_median+1,
                            ifelse(TP_median == 0, 1/(FN_median+1), 
                                   ifelse(FN_median == 0.01, TP_median+1,
                                          ifelse(TP_median == 0.01, 1/(FN_median+1), TP_median / FN_median)))))) %>%
  mutate(
    Disease = paste0(toupper(substr(Disease, 1, 1)), substr(Disease, 2, nchar(Disease))),
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
    # FC_log = log(FC) )
    FC_log = ifelse(FC == 0, NA,
                    ifelse(FC<0, log(-FC), log(FC))),
    neg_log10_p_adj = ifelse(is.na(FC_log), 0, neg_log10_p_adj),
    neg_log10_p_adj1 = ifelse(neg_log10_p_adj > 10, 10, neg_log10_p_adj))

# Correct names
names1 = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/variable_names.csv") %>%
  dplyr::rename(Variable = variable)
pvalue_df_labs1 = pvalue_df_labs1 %>%
  dplyr::left_join(names1) %>%
  dplyr::rename(Variable1 = Variable,
                Variable = output)

# Plot
p <- ggballoonplot(pvalue_df_labs1, x = "Variable", y = "Disease",
                   fill = "FC_log",
                   size = "neg_log10_p_adj1",
                   # size.range = c(1, 10),
                   ggtheme = theme_bw()) +
  scale_size(breaks = c(0, -log10(0.05), 2, 3),
             range = c(1, 7),
             labels = c("ns", "*", "**", "***"),
             limits = c(-log10(0.05), max(pvalue_df_labs1$neg_log10_p_adj1, na.rm = TRUE))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0) +
  guides(size = guide_legend(title="Log10 AdjP", nrow = 3, title.vjust = 0.5),
         fill = guide_colorbar(title="Log10 TP/FN", title.vjust = 0.75)) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
        axis.text.x = element_text(colour="black", angle=45, hjust=1),
        axis.title.x = element_text(size=12, colour="black", face="bold"),
        axis.text.y = element_text(colour="black"),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(colour="black", face="bold"),
        legend.position = "right",
        legend.box="vertical",
        legend.margin=margin()); p
ggsave(plot = p, filename = paste0(results, "/Balloonplot_explainability.png"), width = 8, height = (round(length(unique(pvalue_df_labs1$Disease))/5)+2.5), units = "in", dpi = 300)

